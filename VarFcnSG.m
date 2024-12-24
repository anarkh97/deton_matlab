% The Stiffened gas equation of state. Implementation is based on M2C. Here we port
% the C++ code to Matlab.

classdef VarFcnSG < VarFcnBase
   
%     member variables -- External code can get them but not SET them.
    properties (SetAccess = private)
        gamma double
        pstiff double
%         rho0 double
    end
%     calculated properties.
    properties (SetAccess = private, Hidden)
        invgamma double % 1/gamma
        gamma1 double % gamma - 1
        invgamma1 double % 1/(gamma - 1)
    end
%     public member functions
    methods (Access = public)
%         constructor
        function obj = VarFcnSG(gamma, pstiff)
            arguments
                gamma double
                pstiff double
%                 rho0 double
            end
            
            if nargin ~= 2
                error("*** Error: Required number for inputs for Stiffened gas " ...
                    + "equation of state were not provided.");
            end
            
            obj.gamma = gamma;
            obj.pstiff = pstiff;
%             obj.rho0 = rho0;
            
            obj.invgamma = 1/gamma;
            obj.gamma1 = gamma - 1;
            obj.invgamma1 = 1/obj.gamma1;
            
        end
    end
    
    methods (Access = public) % set access based on VarFcnBase
        function p = GetPressure(obj, rho, e)
            arguments
                obj
                rho double
                e double
            end
            p = obj.gamma1*rho*e - obj.gamma*obj.pstiff;
        end
        function e = GetInternalEnergyPerUnitMass(obj, rho, p)
            arguments
                obj
                rho double
                p double
            end
            e = (p + obj.gamma*obj.pstiff)/(obj.gamma1*rho);
        end
        function rho = GetDensity(obj, p, e)
            arguments
                obj
                p double
                e double
            end
            rho = (p + obj.gamma*obj.pstiff)/(obj.gamma1*e);
        end
        function DpDrho = GetDpDrho(obj, rho, e)
            arguments
               obj
               rho double
               e double
            end
            DpDrho = obj.gamma1*e;
        end
%         here rho and e are unused, added regardless so that we are consistent.
        function Gamma = GetBigGamma(obj, rho, e)
            arguments
               obj
               rho double
               e double
            end
            Gamma = obj.gamma1;
        end
    end
end