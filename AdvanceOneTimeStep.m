function V = AdvanceOneTimeStep(V, varfcn, x, Id, dt)

arguments
   V(3, :) double {mustBeNumeric, mustBeReal}
   varfcn(1, :) VarFcnBase
   x(1, :) double {mustBeNumeric, mustBeReal}
   Id(1, :) uint16 {mustBeNumeric}
   dt double {mustBeReal}
end

if numel(Id) ~= numel(x)
   error("*** Error: Size of arrays Id and x should be identical."); 
end

if numel(varfcn) ~= max(Id)
   error("*** Error: Enough fluid materials were not defined. Number of materials is " ...
       + "%d, whereas the largest material id in the fluid domain is %d", ...
       numel(varfcn), max(Id));
end

N = length(x);
num_of_interfaces = length(varfcn)-1;
dx = 1/N;

% Compute conservative variables
U = V;
for i=1:N
    U(:, i) = varfcn(Id(i)).PrimitiveToConservative(V(:, i));
end

% Evaluate flux term
F = zeros(3, N);
nu = zeros(1, N);

for i=1:N
       
    rho = V(1, i);
    u   = V(2, i);
    p   = V(3, i);
    
    E   = U(3, i);
    
    e   = varfcn(Id(i)).GetInternalEnergyPerUnitMass(rho, p);
%     try
    c   = varfcn(Id(i)).GetSoundSpeed(rho, e);
%     catch ME
%         error("*** Error: Negative speed of sound encountered at Point[%d] x = %e. " ...
%             + "The state variables are (%e, %e, %e) with material id %d.", ...
%             i, x(i), rho, u, p, Id(i));
%     end
    
    F(1, i) = rho*u;
    F(2, i) = rho*u*u + p;
    F(3, i) = u*(E + p);
    
    nu(i)   = 0.5*(abs(u) + c);
    
end
  
nu_interface = max(nu(1:end-1), nu(2:end));
    
% Evaluate numerical flux
F_interm = 0.5*(F(:, 1:end-1) + F(:, 2:end)) + nu_interface.*(U(:, 1:end-1) - U(:, 2:end));

% Calculate flux at material interface
k = 1;
F_intermjm = zeros(3, num_of_interfaces);
F_intermjp = zeros(3, num_of_interfaces);
interface_index = find(diff(Id));
for i = 1:N
    
    if ismember(i, interface_index)
        
        Um = U(:, i);
        Up = U(:, i+1);
        
        Vm = V(:, i);
        Vp = V(:, i+1);
    
%         NOTE: This part of the code only works for Stiffened gas equation of state.
%         TODO: Expand to incorporate other equations as well.
%         Left state
        p0l    = 0; 
        rhol   = Vm(1); 
        ul     = Vm(2); 
        pl     = Vm(3);
        el     = varfcn(Id(i)).GetInternalEnergyPerUnitMass(rhol, pl);
        gammal = varfcn(Id(i)).GetBigGamma(rhol, el) + 1; 
        cl     = varfcn(Id(i)).GetSoundSpeed(rhol, el);
%         Right state
        p0r    = 0; 
        rhor   = Vp(1); 
        ur     = Vp(2); 
        pr     = Vp(3);
        er     = varfcn(Id(i+1)).GetInternalEnergyPerUnitMass(rhor, pr);
        gammar = varfcn(Id(i+1)).GetBigGamma(rhor, er) + 1;
        cr     = varfcn(Id(i+1)).GetSoundSpeed(rhor, er);

        if pl<0 || pr<0 || rhol<0 || rhor<0
            error("*** Error: Negative pressure or density at the material interface %d", ...
                Id(i));
        end
        
%         EXACT RIEMANN SOLVER
        rsolve = rsol_func_mex(gammal, p0l, gammar, p0r, rhol, ul, pl, rhor, ur, pr);
        
%         Left star states 
        psl = rsolve(2);
        usl = rsolve(4);  
        rhosl = rsolve(6);
%         Right star states
        psr = rsolve(1); 
        usr = rsolve(3);
        rhosr = rsolve(5); 
        
        esl = varfcn(Id(i)).GetInternalEnergyPerUnitMass(rhosl, psl);
        esr = varfcn(Id(i+1)).GetInternalEnergyPerUnitMass(rhosr, psr);
        
        Esl = rhosl*(esl + 0.5*usl*usl);
        Esr = rhosr*(esr + 0.5*usr*usr);

%         Left conservative state
        Usl = [rhosl; rhosl*usl; Esl];
%         Right conservative state
        Usr = [rhosr; rhosr*usr; Esr];
        
%         Left flux term
        Fsl = [rhosl*usl; rhosl*usl*usl + psl; usl*(Esl + psl)];
%         Right flux term
        Fsr = [rhosr*usr; rhosr*usr*usr + psr; usr*(Esr + psr)];
        
        csl = varfcn(Id(i)).GetSoundSpeed(rhosl, esl);
        csr = varfcn(Id(i+1)).GetSoundSpeed(rhosr, esr);
        
        nusl = 0.5*max(abs(ul) + cl, abs(usl) + csl);
        nusr = 0.5*max(abs(ur) + cr, abs(usr) + csr);

        F_intermjm(:, k) = 0.5*(F(:, i) + Fsl) + nusl*(Um - Usl);
        F_intermjp(:, k) = 0.5*(Fsr + F(:, i+1)) + nusr*(Usr - Up);
        
        k = k + 1;
    end
    
end

% Append boundary values   
F = [F(:, 1) , F_interm , F(:, end)];

% Forward Euler time stepping
k = 1;
% dFDebug = zeros(3, N);
for i = 1:N
    dF = zeros(3, 1);
    if ismember(i, interface_index)
%         Left side of the interface
        dF(:) = (F(:, i) - F_intermjm(:, k))/dx;
    elseif ismember(i, interface_index+1)
%         Right side of the interface
        dF(:) = (F_intermjp(:, k) - F(:, i+1))/dx;
        k = k + 1;
    else
%         Not near the interface
        dF(:) = (F(:, i) - F(:, i+1))/dx;
    end
%     dFDebug(:, i) = dF(:);
    U(:, i) = U(:, i) + dt*dF;
    
    V(:, i) = varfcn(Id(i)).ConservativeToPrimitive(U(:, i));
end
end