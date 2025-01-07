function V = AdvanceOneTimeStep(V, mixture, x, dt)

arguments
   V(5, :) double {mustBeNumeric, mustBeReal}
   mixture MixtureModel
   x(1, :) double {mustBeNumeric, mustBeReal}
   dt double {mustBeReal}
end

if numel(x) ~= numel(V(1, :))
   error("*** Error: Size of mesh does not match with the size of state variables."); 
end

N = length(x);
dx = 1/N;

% Compute conservative variables
U = V;
for i=1:N
    U(:, i) = mixture.PrimitiveToConservative(V(:, i));
end

% Evaluate flux term
F = zeros(5, N);
nu = zeros(1, N);

for i=1:N
       
    rho_1  = V(1, i);
    rho_2  = V(2, i);
    u      = V(3, i);
    p      = V(4, i);
    lambda = V(5, i);
    
    E   = U(4, i);
    
    rho = lambda*rho_1 + (1 - lambda)*rho_2;
    e   = mixture.GetInternalEnergyPerUnitMass(rho_1, rho_2, p, lambda);
    try
        c   = mixture.GetSoundSpeed(rho_1, rho_2, e, lambda);
    catch ME
        fprintf("%s\n", ME.message);
        error("*** Error: Negative speed of sound encountered at Point[%d] x = %e. " ...
            + "The state variables are (%e, %e, %e).", ...
            i, x(i), rho, u, p);
    end
    
    F(1, i) = lambda*rho_1*u;
    F(2, i) = (1 - lambda)*rho_2*u;
    F(3, i) = rho*u*u + p;
    F(4, i) = u*(E + p);
    F(5, i) = lambda*u;
    
    nu(i)   = abs(u) + c;
    
end
    
% Evaluate numerical flux -- LLF
nu = max(nu(1:end-1), nu(2:end));
F_interm = 0.5*(F(:, 1:end-1) + F(:, 2:end) - nu.*(U(:, 2:end) - U(:, 1:end-1)));

% Append boundary values   
F = [F(:, 1) , F_interm , F(:, end)];

% Compute source terms (lambda advection)
S = zeros(5, N);
u_x = [V(3, 1), V(3, :), V(3, end)];
S(5, :) = (u_x(3:end) - u_x(1:end-2))/(2*dx);

% Forward Euler time stepping
for i = 1:N
    dF = (F(:, i) - F(:, i+1))/dx;
%     dFDebug(:, i) = dF(:);
    U(:, i) = U(:, i) + dt*dF + dt*S(:, i);
    
    V(:, i) = mixture.ConservativeToPrimitive(U(:, i));
end
end