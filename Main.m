clear; 
close all;
clc;

% This code simulates the SG-VDW example case from https://doi.org/10.1006/jcph.2002.7143 
% Section 8.3. We use the five equation model, which tracks the fractional density of 
% each fluid material. The volume fraction of fluid material is advected in space and time to
% calculate relevant pressure and internal energy in the control volume.
% NOTE: This code is intended to work only for two fluid materials. The code will
% be expanded further to incorporate a reaction equation, and eventually simulate 
% 1-D detonations.

% units: kg, m, s, Pa

% Mixture model -- has two fluid materials (i.e. VarFcn's).
mixture = MixtureModel(VarFcnSG(4.4, 6e8), VarFcnVDW(1.4, 5, 1e-3, 0));

% Initial Conditions
% User inputs
rhoIC = [1000 50];  % density 
uIC = [0 0];        % velocity
pIC = [1e9 1e5];    % pressure
lambdaIC = [1-eps eps];   % volume fraction

% Spatial domain
N = 50;
x_in = 0; x_out = 1;
dx = 1/N;
x = (x_in + dx/2):dx:(x_out-dx/2);

% Temporal parameters
t      = 0; 
dt     = 1e-5; 
tf     = 1.0e-2;
n      = 0;
nsteps = round(tf/dt);

% Interface location
% `x_int` stores the x co-ordinate of the material interface.
x_interface = 0.7;

start_time = tic;

% Setup initial state matrix
% Only 5 primitive variables; [\rho_1, \rho_2, u_x, p, \lambda]
Vn = zeros(5, N);
% Not used. Kept here for reference.
% Only 5 conservative variables; [\rho_1*\lambda, \rho_2*(1-\lambda), \rho u_x, E]
% Un = zeros(5, N);

% Populate primitive variables and material ids
for i=1:N
   
    if x(i) <= x_interface
%        Left of interface
        Vn(1, i)   = rhoIC(1);
        Vn(2, i)   = rhoIC(2); %eps;
        Vn(3, i)   = uIC(1);
        Vn(4, i)   = pIC(1);
        Vn(5, i)   = lambdaIC(1); % + sign(x(i) - x_interface)*eps;
    else
%        Right of interface
        Vn(1, i)   = rhoIC(1); %eps;
        Vn(2, i)   = rhoIC(2);
        Vn(3, i)   = uIC(2);
        Vn(4, i)   = pIC(2);
        Vn(5, i)   = lambdaIC(2); % + sign(x(i) - x_interface)*eps;
    end
    
end

% For post-processing
% indices for accessing state variables from post processing matrices
keys = reshape(1:5*N, 5, N);
VVn = zeros(5*N, nsteps+1);
VVn(:, 1) = reshape(Vn, [5*N, 1]);

% =================================================================================================
% Time loop
% =================================================================================================
while t<=tf

    n = n + 1;
    t = t + dt;
    
    fprintf("Step %d: t = %e, dt = %e. Computation time: %.4e s.\n", n, t, dt, toc(start_time));
        
%     Advance in time (uses forward euler method and Exact Riemann solver)
    Vn  = AdvanceOneTimeStep(Vn, mixture, x, dt);
    
%     Save data and update loop
    VVn(:, n) = reshape(Vn, [5*N, 1]);

end

% Post processing (t = tf)
close all

rhot  = VVn(keys(1, :), :);
ut    = VVn(keys(2, :), :);
pt    = VVn(keys(3, :), :);
t_vec = linspace(0, t, nsteps+1);

% Density plot
figure(1)
grid on;
box on;
hold on;
plot(x, rhot(:, end), 'LineWidth', 2);
% for k = 1:num_of_interfaces
%     xline(XXn(k, end),'--', 'linewidth', 1)
% end
hold off;
title('Density')
txt = ['t = ' num2str(t_vec(end))];
text(0.1,0.1,txt, 'Units', 'normalized');

% Pressure plot
figure(2)
grid on;
box on;
hold on;
plot(x, pt(:, end), 'LineWidth', 2);
% for k = 1:num_of_interfaces
%     xline(XXn(k, end),'--', 'linewidth', 1)
% end
hold off;
title('Pressure')
txt = ['t = ' num2str(t_vec(end))];
text(0.1,0.1,txt, 'Units', 'normalized');

% Velocity plot
figure(3)
grid on;
box on;
hold on;
plot(x, ut(:, end), 'LineWidth', 2);
% for k = 1:num_of_interfaces
%     xline(XXn(k, end),'--', 'linewidth', 1)
% end
hold off;
title('Velocity')
txt = ['t = ' num2str(t_vec(end))];
text(0.1,0.1,txt, 'Units', 'normalized');

% Post processing (animation)
close all
figure(4)
ff = 20; % frames per second
for n = 1:ff:nsteps
% % % Toggle comments to plot different variables % % %
%     plot(x, rhot(:, n), 'linewidth', 2); a = 1;
%     plot(x, ut(:, n), 'linewidth', 2); a = 2;
    plot(x, pt(:, n), 'linewidth', 2); a = 3;
%     for k = 1:num_of_interfaces
%         xline(XXn(k, n),'--', 'linewidth', 1)
%     end
    txt = ['t = ' num2str(t_vec(n))];
    text(0.1,0.1,txt, 'Units', 'normalized');
    
    if a==1
        title('Density')
    elseif a==2
        title('Velocity')
    elseif a==3
        title ('Pressure')
    end
    drawnow
end