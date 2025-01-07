% The extended version of Mie-Gr\:uneisen equation of state. 
% Implementation is based on M2C. Here we port the C++ code to Matlab.
% We are gonna skip temperature computation for now.

classdef VarFcnMGExt < VarFcnBase
   
%     member variables -- External code can get them but not SET them.
    properties (SetAccess = private)
        rho0 double
        c0 double
        Gamma0 double
        s double 
        e0 double 
    end
%     calculated properties.
    properties (SetAccess = private, Hidden)
        eta_min double
        p_min double
        eta_max double
        
        rho0_c0_c0 double % rho0*c0*c0
        Gamma0_over_2 double % Gamma0/2
        Gamma0_rho0 double % Gamma*rho0
        er_min double
        
%         h0 double % reference enthalpy at T0
    end
    methods (Static)
       
        function name = GetEOSName()
            name = "MG Ext. EOS";
        end
         
    end
%     public member functions
    methods (Access = public)
%         constructor
        function obj = VarFcnMGExt(rho0, c0, Gamma0, s, e0, eta_min)
            arguments
                rho0 double
                c0 double
                Gamma0 double
                s double
                e0 double
%                 h0 double
                eta_min double
            end
            
            if nargin ~= 7
                error("*** Error: Required number for inputs for extended Mie-Gruneisen equation of state were not provided.");
            end
            
            obj.rho0 = rho0;
            obj.c0 = c0;
            obj.Gamma0 = Gamma0;
            obj.s = s;
            obj.e0 = e0;
%             obj.h0 = h0;
            obj.eta_min = eta_min;
            
            obj.p_min = rho0*c0*c0*eta_min;
            obj.eta_max = 1/s - 1e-6;
            
            obj.rho0_c0_c0 = rho0*c0*c0;
            obj.Gamma0_over_2 = 0.5*Gamma0;
            obj.Gamma0_rho0 = Gamma0*rho0;
            obj.er_min = 0.5*c0*c0*eta_min*eta_min + e0;
                        
        end
    end
    
    methods (Access = public) % set access based on VarFcnBase
        function p = GetPressure(obj, rho, e)
            arguments
                obj
                rho double
                e double
            end
            eta = 1.0 - obj.rho0/rho;
            if eta>=0.0
                p = obj.rho0_c0_c0*eta*(1.0 - obj.Gamma0_over_2*eta)/((1.0-obj.s*eta)*(1.0-obj.s*eta)) ...
                    + obj.Gamma0_rho0*(e-obj.e0);
                return;
            elseif eta>=obj.eta_min
                p = obj.rho0_c0_c0*eta*(1.0 - obj.Gamma0_over_2*eta) + obj.Gamma0_rho0*(e-obj.e0);
                return;
            else
                p = obj.rho0_c0_c0*obj.eta_min*(1.0 - obj.Gamma0_over_2*(2.0*eta-obj.eta_min)) ...
                    + obj.Gamma0_rho0*(e-obj.e0);
                return;
            end
        end
        function e = GetInternalEnergyPerUnitMass(obj, rho, p)
            arguments
                obj
                rho double
                p double
            end
            eta = 1.0 - obj.rho0/rho;
            e = (p - GetPr(eta))/obj.Gamma0_rho0 + GetEr(eta); 
        end
        function rho = GetDensity(obj, p, e)
            arguments
                obj
                p double
                e double
            end
            
            ptrial = obj.GetPressure(obj.rho0, e);
            if ptrial == p
                rho = obj.rho0;
                return;
            end
            if ptrial<p 
                c = p - obj.Gamma0_rho0*(e - obj.e0);
                a = c*obj.s*obj.s + obj.Gamma0_over_2*obj.rho0_c0_c0;
                b = -(2.0*c*obj.s + obj.rho0_c0_c0);
                
                if(a==0) % linear equation: eta = -c/b, rho = rho0/(1-eta)
                    rho = obj.rho0/(1 + c/b);
                    return;
                end
                
                b2m4ac = b*b - 4.0*a*c;
                if (b2m4ac<0) 
                    error("*** Error: The Extended M-G EOS is invalid for the given p(%e) and e(%e) " ...
                        + "--- unable to solve it for rho.", p, e);
                end

                b2m4ac = sqrt(b2m4ac);
                eta1 = (-b + b2m4ac)/(2.0*a);
                eta2 = (-b - b2m4ac)/(2.0*a);
                
                if(eta1 > eta2)
                    % eta2 should be the bigger one
                    [eta1, eta2] = deal(eta1, eta2);
                end
                
                if eta1>=1.0 || eta2<0  
                    % both are negative, or both are positive
                    error("*** Error: Cannot determine the solution (rho) of the Extended M-G EOS " ...
                        + "(eta1 = %e, eta2 = %e).", eta1, eta2);
                end
                
                if eta2 >= 1.0
                    % eta1 should be the valid root.
                    if eta1 < 0.0
                       error("*** Error: Cannot find the solution (rho) of the Extended M-G EOS " ...
                           + "(eta1 = %e, eta2 = %e).", eta1, eta2);
                    end
                    rho = obj.rho0/(1.0 - eta1);
                    return;
                end
                
                % now, eta2 must be between 0 and 1 (i.e., valid)
                if eta1>=0.0 && eta1 ~= eta2
                    error("*** Error: Cannot find the solution (rho) of the Extended M-G EOS " ...
                        + "(eta1 = %e, eta2 = %e).", eta1, eta2);
                end
                
                rho = obj.rho0/(1.0 - eta2);
                return;
            end
            
            % If we get here, it means ptrial>p ==> eta<0
            c = (p - obj.Gamma0_rho0*(e-obj.e0))/obj.rho0_c0_c0;
            eta = (1.0-c/obj.eta_min)/obj.Gamma0 + 0.5*obj.eta_min;
            
            if eta <= obj.eta_min
               rho = obj.rho0/(1 - eta);
               return;
            end
            
            % If we get here, it means eta_min < eta < 0
            
            %solving a quadratic equation a*eta^2 + b*eta + c = 0 for eta ==> rho = rho0/(1-eta)
            b2m4ac = 1.0 - 2.0*obj.Gamma0*c;
            
            if b2m4ac<0 
                error("*** Error: The Extended M-G EOS is invalid for the given p(%e) and e(%e) " ...
                    + "--- unable to solve it for rho (eta_min<eta<0).", p, e);
            end
            
            b2m4ac = sqrt(b2m4ac);
            eta = (1.0 - sqrt(b2m4ac))/obj.Gamma0;

            if eta>0 || eta<obj.eta_min 
               error("*** Error: Cannot determine the solution (rho) of the Extended M-G EOS " ...
                   + "(eta = %e).", eta);
            end

            rho = obj.rho0/(1.0 - eta);
  
        end
        function Pi = GetBigPi(obj, rho)
            arguments
               obj
               rho double
            end
            error("*** Error: VarFcnMGExt::GetBigPi() not implemented.\n");
        end
%         e is unused but added to be consistent
        function Gamma = GetBigGamma(obj, rho, e)
            arguments
               obj
               rho double
               e double
            end
            Gamma = obj.Gamma0_rho0/rho;
        end
        function DpDrho = GetDpDrho(obj, rho, e)
            arguments
                obj
                rho double
                e double
            end
            eta = 1 - obj.rho0/rho;
            DpDrho = (obj.GetDPrDeta(eta) - obj.Gamma0_rho0*obj.GetDErDeta(eta))*obj.GetDetaDrho(rho);
        end
    end
%     private member functions
    methods (Access = private)     
        function DPrDeta = GetDPrDeta(obj, eta)
           arguments
               obj
               eta double
           end
           if eta >= 0.0 
               if eta >= obj.eta_max
                  fprintf("*** Warning: (Ext-MG EOS) Exceeding the compression limit.");
                  eta = obj.eta_max;
               end
               f = 1.0/(1.0-obj.s*eta);
               DPrDeta = obj.rho0_c0_c0*(1+obj.s*eta)*f*f*f;
               return;
           elseif eta >= obj.eta_min
               DPrDeta = obj.rho0_c0_c0;
               return;
           else
               DPrDeta = 0.0;
               return;
           end
        end
        function DErDeta = GetDErDeta(obj, eta)
           arguments
               obj
               eta double
           end
           if eta >= 0.0
               if eta >= obj.eta_max
                  fprintf("*** Warning: (Ext-MG EOS) Exceeding the compression limit.");
                  eta = obj.eta_max;
               end
               f = 1.0/(1.0-obj.s*eta);
               DErDeta = obj.c0*obj.c0*eta*f*f*f;
               return;
           elseif eta >= obj.eta_min
               DErDeta = obj.c0*obj.c0*eta;
               return;
           else
               DErDeta = obj.p_min/obj.rho0;
               return;
           end
        end
        function DetaDrho = GetDetaDrho(obj, rho)
           arguments
               obj
               rho double
           end
           DetaDrho = obj.rho0/(rho*rho);
        end
        function pr = GetPr(obj, eta)
           arguments
               obj
               eta double
           end
           if eta >= 0.0 
               if eta >= obj.eta_max
                  fprintf("*** Warning: (Ext-MG EOS) Exceeding the compression limit.");
                  eta = obj.eta_max;
               end
               f = 1.0/(1.0-obj.s*eta);
               pr = obj.rho0_c0_c0*eta*f*f;
               return;
           elseif eta >= obj.eta_min
               pr = obj.rho0_c0_c0*eta;
               return;
           else
               pr = obj.p_min;
               return;
           end
        end
        function er = GetEr(obj, eta)
           arguments
               obj
               eta double
           end
           if eta >= 0.0
               if eta >= obj.eta_max
                  fprintf("*** Warning: (Ext-MG EOS) Exceeding the compression limit.");
                  eta = obj.eta_max;
               end
               f = 1.0/(1.0-obj.s*eta);
               er = 0.5*obj.c0*obj.c0*eta*eta*f*f + obj.e0;
               return;
           elseif eta >= obj.eta_min
               er = 0.5*obj.c0*obj.c0*eta*eta + obj.e0;
               return;
           else
               er = obj.er_min + obj.p_min/obj.rho0*(eta - obj.eta_min);
               return;
           end
        end
    end
end