% Base class for equations of state. Implementation is based on M2C. Here we port
% the C++ code to Matlab. We ignore some functions and obly implement the ones we need.

classdef (Abstract) VarFcnBase < matlab.mixin.Heterogeneous
    
    methods (Access = public)
%         State conversion methods
        function U = PrimitiveToConservative(obj, V)
            % Performing a primitive to conservative operation
            arguments
                obj
%                 mustBeNumeric might be redundant
                V(3, 1) double {mustBeNumeric, mustBeReal}
            end

            % Primitive state variables
            % V[1] -> \rho
            % V[2] -> u_x
            % V[3] -> p

            % Conservative state variables
            % U[1] -> \rho
            % U[2] -> \rho u_x
            % U[3] -> E

            U = V;

            U(1) = V(1);
            U(2) = V(1)*V(2);
            
            e = obj.GetInternalEnergyPerUnitMass(V(1), V(3)); % pass rho and p
            U(3) = V(1)*(e + 0.5*V(2)*V(2));
        end
        function V = ConservativeToPrimitive(obj, U)
            % Performing a conservative to primitive operation
            arguments
                obj
%                 mustBeNumeric might be redundant
                U(3,1) double {mustBeNumeric, mustBeReal}
            end

            % Primitive state variables
            % V[1] -> \rho
            % V[2] -> u_x
            % V[3] -> p

            % Conservative state variables
            % U[1] -> \rho
            % U[2] -> \rho u_x
            % U[3] -> E

            V = U;
            invrho = 1/U(1);
            V(1) = U(1);
            V(2) = U(2)*invrho;
            
            e = (U(3) - 0.5*V(1)*V(2)*V(2))*invrho;
            V(3) = obj.GetPressure(V(1), e); % pass rho and e
        end
        function c = GetSoundSpeed(obj, rho, e)
            arguments
               obj
               rho double
               e double
            end
            c2 = obj.GetDpDrho(rho, e) + obj.GetPressure(rho, e)/rho*obj.GetBigGamma(rho, e);
            if c2 <= 0
                error("*** Error: Cannot calculate speed of sound (Square-root of a negative number)")
            end
            c = sqrt(c2);
        end
    end
    
%     Abstract functions
    methods (Abstract, Access = public)
        p = GetPressure(obj, rho, e)
        e = GetInternalEnergyPerUnitMass(obj, rho, p)
        rho = GetDensity(obj, p, e)
        Gamma = GetBigGamma(obj, rho, e)
        DpDrho = GetDpDrho(obj, rho, e)
    end
    
end