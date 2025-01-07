% This model represents the mixture equation of state. It takes
% two different equations of state and the volume fraction parameter
% which should present material 1 when 1 and material 2 when 0.

% Mixture properties are given as,
%   1. Density         : \rho = \rho_1*\lambda + \rho_2*(1 - \lambda)
%   2. Internal energy : \rho*e = \rho_1*e_1*\lambda + \rho_2*e_2*(1 - \lambda)
%   3. Pressure        : isobaric i.e., p_1 = p_2 = p

classdef MixtureModel
   
    properties (SetAccess = private)
       
        mat1 VarFcnBase
        mat2 VarFcnBase
        
    end
    
    methods (Access = public)
        function obj = MixtureModel(mat_1, mat_2)
            arguments
                mat_1 VarFcnBase
                mat_2 VarFcnBase
            end
            
            obj.mat1 = mat_1;
            obj.mat2 = mat_2;
            
        end
%         State conversion methods
        function U = PrimitiveToConservative(obj, V)
            % Performing a primitive to conservative operation
            arguments
                obj
%                 mustBeNumeric might be redundant
                V(5, 1) double {mustBeNumeric, mustBeReal}
            end

            % Primitive state variables
            % V[1] -> \rho_1
            % V[2] -> \rho_2
            % V[3] -> u_x
            % V[4] -> p
            % V[5] -> \lambda

            % Conservative state variables
            % U[1] -> \rho_1*\lambda
            % U[2] -> \rho_2*(1 - \lambda)
            % U[3] -> \rho u_x
            % U[4] -> E
            % U[5] -> \lambda

            U = V;

            U(1) = V(5)*V(1);
            U(2) = (1 - V(5))*V(2);
            
            rho  = U(1) + U(2);
            U(3) = rho*V(3);
            
            e = obj.GetInternalEnergyPerUnitMass(V(1), V(2), V(4), V(5)); % pass rho_1, rho_2, p and lambda
            U(4) = rho*(e + 0.5*V(3)*V(3));
            
            U(5) = V(5);
        end
        function V = ConservativeToPrimitive(obj, U)
            % Performing a conservative to primitive operation
            arguments
                obj
%                 mustBeNumeric might be redundant
                U(5, 1) double {mustBeNumeric, mustBeReal}
            end

            % Primitive state variables
            % V[1] -> \rho_1
            % V[2] -> \rho_2
            % V[3] -> u_x
            % V[4] -> p
            % V[5] -> \lambda

            % Conservative state variables
            % U[1] -> \rho_1*\lambda
            % U[2] -> \rho_2*(1 - \lambda)
            % U[3] -> \rho u_x
            % U[4] -> E
            % U[5] -> \lambda

            V = U;
            
            rho = U(1) + U(2);
            invrho = 1/rho;
        
%             Avoid division by zero
%             if U(5) < eps
%                 V(1) = eps;
%             elseif abs(1 - U(5)) < eps
%                 V(2) = eps;
%             else
%                 V(1) = U(1)/U(5);
%                 V(2) = U(2)/(1 - U(5));
%             end
            V(1) = U(1)/U(5);
            V(2) = U(2)/(1 - U(5));
            V(3) = U(3)*invrho;
            
            e = (U(4) - 0.5*rho*V(3)*V(3))*invrho;
            V(4) = obj.GetPressure(V(1), V(2), e, V(5)); % pass rho_1, rho_2, e and lambda
            
            V(5) = U(5);
        end
        function c = GetSoundSpeed(obj, rho_1, rho_2, e, lambda)
            arguments
               obj
               rho_1 double
               rho_2 double
               e double
               lambda double
            end

            rho  = lambda*rho_1 + (1 - lambda)*rho_2;
            p    = obj.GetPressure(rho_1, rho_2, e, lambda);
            e_1  = obj.mat1.GetInternalEnergyPerUnitMass(rho_1, p);
            e_2  = obj.mat2.GetInternalEnergyPerUnitMass(rho_2, p);
            
            try
                c_1 = obj.mat1.GetSoundSpeed(rho_1, e_1);
            catch ME
                fprintf("%s\n", ME.message);
                error("rho_1 = %e, e_1 = %e, p = %e\n", ...
                    rho_1, e_1, p);
            end
            try
                c_2 = obj.mat2.GetSoundSpeed(rho_2, e_2);
            catch ME
                fprintf("%s\n", ME.message);
                error("rho_2 = %e, e_2 = %e, p = %e\n", ...
                    rho_2, e_2, p);
            end

            Gamma1 = obj.mat1.GetBigGamma(rho_1, e_1);
            Gamma2 = obj.mat2.GetBigGamma(rho_2, e_2);
            
%             TODO: Need to check if these relations are valid.
            Xi1    = 1/Gamma1;
            Xi2    = 1/Gamma2;
            Xi     = lambda*Xi1 + (1 - lambda)*Xi2;
            
            y1     = lambda*rho_1/rho;
            y2     = (1 - lambda)*rho_2/rho;
            
            c2     = (y1*Xi1*c_1*c_1 + y2*Xi2*c_2*c_2)/Xi;
            if c2 <= 0
                error("*** Error: Cannot calculate speed of sound for the mixture " ...
                    + "(Square-root of a negative number)")
            end
            c      = sqrt(c2);
            
        end
        
        function p = GetPressure(obj, rho_1, rho_2, e, lambda)
            arguments
                obj
                rho_1 double
                rho_2 double
                e double
                lambda double
            end
            
%             We solve the pressure equation iteratively. Solution should exist as for general
%             equation of state (e.g., Mie_Gr\:uneisen) pressure can be defined explicitly in terms
%             for fractional densities and internal energy (e). 
%             NOTE: With iterative approach we can use existing implementations of VarFcn's w/o
%             adding new functinalities.
%             This might slow down the code.
            try
                if abs(1 - lambda) < 1e-8
                    % found a pure fluid, save computation
                    p = obj.mat1.GetPressure(rho_1, e);
                elseif abs(lambda) < 1e-8
                    % found a pure fluid, save computation
                    p = obj.mat2.GetPressure(rho_2, e);
                else
                    eqn     = @(x) obj.Fun(rho_1, rho_2, e, lambda, x);
                    options = optimoptions('fsolve', 'Display', 'none', 'MaxIterations', 20);
                    p = fsolve(eqn, 0, options);
                end
            catch ME
                error("rho_1: %e, rho_2: %e, e: %e, lambda: %e\n", rho_1, rho_2, e, lambda);
            end
            
        end
        function e = GetInternalEnergyPerUnitMass(obj, rho_1, rho_2, p, lambda)
            arguments
                obj
                rho_1 double
                rho_2 double
                p double
                lambda double
            end
            
            e_1 = obj.mat1.GetInternalEnergyPerUnitMass(rho_1, p);
            e_2 = obj.mat2.GetInternalEnergyPerUnitMass(rho_2, p);
            
            rho    = lambda*rho_1 + (1 - lambda)*rho_2;
            e      = (lambda*rho_1*e_1 + (1 - lambda)*rho_2*e_2)/rho;
            
        end
        
    end
    
    methods (Access = private)
%        \rho*e = \rho_1*e_1*\lambda + \rho_2*e_2*(1 - \lambda)
        function r = Fun(obj, rho_1, rho_2, e, lambda, p)
            arguments
                obj
                rho_1 double
                rho_2 double
                e double
                lambda double
                p double
            end
            
            rho    = lambda*rho_1 + (1 - lambda)*rho_2;
            
            term1  = lambda*rho_1*obj.mat1.GetInternalEnergyPerUnitMass(rho_1, p);
            term2  = (1 - lambda)*rho_2*obj.mat2.GetInternalEnergyPerUnitMass(rho_2, p);
            
            r = rho*e - term1 - term2;
        end
        
    end
    
end