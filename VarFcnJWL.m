% The Jones-Wilkins-Lee equation of state. Implementation is based on M2C. Here we port
% the C++ code to Matlab.

classdef VarFcnJWL < VarFcnBase
   
%     member variables -- External code can get them but not SET them.
    properties (SetAccess = private)
        omega double
        A1 double
        A2 double
        R1 double 
        R2 double 
        rho0 double
    end
%     calculated properties
    properties (SetAccess = private, Hidden)
        R1rho0 double % R1*rho0
        R2rho0 double % R2*rho0
        omega_over_R1rho0 double % omega/(R1/rho0)
        omega_over_R2rho0 double % omega/(R2/rho0)
    end
%     public member functions
    methods (Access = public) 
%         constructor
        function obj = VarFcnJWL(omega, A1, A2, R1, R2, rho0)
            arguments
                omega double
                A1 double
                A2 double
                R1 double
                R2 double
                rho0 double
            end
            
            if nargin ~= 6
                error("*** Error: Required number for inputs for JWL equation of state were not provided.");
            end
            
            obj.omega = omega;
            obj.A1 = A1;
            obj.A2 = A2;
            obj.R1 = R1;
            obj.R2 = R2;
            obj.rho0 = rho0;
            
            obj.R1rho0 = R1*rho0;
            obj.R2rho0 = R2*rho0;
            obj.omega_over_R1rho0 = omega/obj.R1rho0;
            obj.omega_over_R2rho0 = omega/obj.R2rho0;
            
        end
    end
    methods (Access = public) % set access based on VarFcnBase
        function p = GetPressure(obj, rho, e)
            arguments
                obj
                rho double
                e double
            end
            p = obj.omega*rho*e + obj.Fun(rho);
        end
        function e = GetInternalEnergyPerUnitMass(obj, rho, p)
            arguments
                obj
                rho double
                p double
            end
            e = (p-obj.Fun(rho))/(obj.omega*rho);
        end
        function rho = GetDensity(obj, p, e)
            arguments
                obj
                p double
                e double
            end
%             Non-linear equation in density; p - \omega\rho e - f(\rho) = 0
%             where, f() is the non-linear density term in JWL equation of state.
            DensityEquation = @(rho_) p ...
                                     - obj.omega*e ...
                                     - obj.A1*(1.0-obj.omega_over_R1rho0*rho_)*exp(-obj.R1rho0/rho_) ...
                                     - obj.A2*(1.0-obj.omega_over_R2rho0*rho_)*exp(-obj.R2rho0/rho_);
            rho = fsolve(DensityEquation, obj.rho0);
        end
        function DpDrho = GetDpDrho(obj, rho, e)
            arguments
               obj
               rho double
               e double
            end
            DpDrho = obj.omega*e ...
                + obj.A1*(-obj.omega_over_R1rho0 + obj.R1rho0/(rho*rho) - obj.omega_over_R1rho0*obj.R1rho0/rho)*exp(-obj.R1rho0/rho) ...
                + obj.A2*(-obj.omega_over_R2rho0 + obj.R2rho0/(rho*rho) - obj.omega_over_R2rho0*obj.R2rho0/rho)*exp(-obj.R2rho0/rho);
        end
%         here rho and e are unused, added regardless so that we are consistent.
        function Gamma = GetBigGamma(obj, rho, e)
            arguments
               obj
               rho double
               e double
            end
            Gamma = obj.omega;
        end
    end
%     private member functions
    methods (Access = private)
        function ret = Fun(obj, rho)
           arguments
               obj
               rho double
           end
           ret = obj.A1*(1-obj.omega_over_R1rho0*rho)*exp(-obj.R1rho0/rho) ...
               + obj.A2*(1-obj.omega_over_R2rho0*rho)*exp(-obj.R2rho0/rho);
        end
    end
end