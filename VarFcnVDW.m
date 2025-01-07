% The van der Waals equation of state. Implementation is based on 10.1002/ﬂd

classdef VarFcnVDW < VarFcnBase
   
%     member variables -- External code can get them but not SET them.
    properties (SetAccess = private)
        gamma double
        a double
        b double
        c double
    end
%     calculated properties.
    properties (SetAccess = private, Hidden)
        invgamma double % 1/gamma
        gamma1 double % gamma - 1
    end
    methods (Static)
       
        function name = GetEOSName()
            name = "VDW EOS";
        end
         
    end
%     public member functions
    methods (Access = public)
%         constructor
        function obj = VarFcnVDW(gamma, a, b, c)
            arguments
                gamma double
                a double
                b double
                c double
            end
            
            if nargin ~= 4
                error("*** Error: Required number for inputs for van der Waals " ...
                    + "equation of state were not provided.");
            end
            
            obj.gamma = gamma;
            obj.a = a;
            obj.b = b;
            obj.c = c;
            
            obj.invgamma = 1/gamma;
            obj.gamma1 = gamma - 1;
            
        end
    end
    
    methods (Access = public) % set access based on VarFcnBase
        function p = GetPressure(obj, rho, e)
            arguments
                obj
                rho double
                e double
            end
            factor = obj.gamma1/(1 - obj.b*rho);
            p = factor*(rho*e + obj.a*rho*rho - obj.c) - (obj.a*rho*rho + obj.c);
        end
        function e = GetInternalEnergyPerUnitMass(obj, rho, p)
            arguments
                obj
                rho double
                p double
            end
            factor = obj.gamma1/(1 - obj.b*rho);
            p_tot = p + obj.a*rho*rho + obj.c;
            rho_tot = rho*factor;
            e = p_tot/rho_tot + obj.c/rho - obj.a*rho;
        end
        function rho = GetDensity(obj, p, e)
            arguments
                obj
                p double
                e double
            end
            error("*** Error: VarFcnVDW::GetDensity() has not been implemented.");
        end
        function DpDrho = GetDpDrho(obj, rho, e)
            arguments
               obj
               rho double
               e double
            end
            factor = obj.gamma1/(1 - obj.b*rho);
            factor1 = obj.b*factor/(1 - obj.b*rho);
            term1 = factor1*(rho*e + obj.a*rho*rho - obj.c);
            term2 = factor*(e + 2*obj.a*rho);
            DpDrho = term1 + term2 - 2*obj.a*rho;
        end
%         here e is unused, added regardless so that we are consistent.
        function Gamma = GetBigGamma(obj, rho, e)
            arguments
               obj
               rho double
               e double
            end
            Gamma = obj.gamma1/(1 - obj.b*rho);
        end
        function Pi = GetBigPi(obj, rho)
            arguments
               obj
               rho double
            end
            
            factor = obj.gamma1/(1 - obj.b*rho);
            Pi = (1 - factor)*obj.a*rho*rho + (factor + 1)*obj.c;
        end
    end
end