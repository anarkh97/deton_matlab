function [dt, cfl] = CalculateCFLAndStepSize(dx, dt, cflmax, mixture, V)

arguments
   dx double {mustBeNumeric, mustBeReal}
   dt double {mustBeNumeric, mustBeReal}
   cflmax double {mustBeNumeric, mustBeReal}
   mixture MixtureModel
   V(5, :) double {mustBeNumeric, mustBeReal}
end

N = length(V);

% compute characteristic speed
s = -realmax;
for i=1:N
       
    rho_1  = V(1, i);
    rho_2  = V(2, i);
    u      = V(3, i);
    p      = V(4, i);
    lambda = V(5, i);
    
    rho = lambda*rho_1 + (1 - lambda)*rho_2;
    e   = mixture.GetInternalEnergyPerUnitMass(rho_1, rho_2, p, lambda);
    try
        c   = mixture.GetSoundSpeed(rho_1, rho_2, e, lambda);
    catch ME
        fprintf("%s\n", ME.message);
        error("*** Error: Negative speed of sound. The state variables are (%e, %e, %e).", ...
            rho, u, p);
    end
    
    s1 = abs(u - c);
    s2 = abs(u);
    s3 = abs(u + c);
    
    s_temp = max([s1, s2, s3]);
    
    if(s_temp >= s)
        s = s_temp;
    end
    
end

% calculate cfl number
cfl = s*dt/dx;

% calculate step size
if cfl > cflmax
   
    cfl = cflmax;
    dt = cflmax*dx/s;
    
end

end

