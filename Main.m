clear; 
close all;
clc;

% This code simulates the PG-PG example case from Arthur Rallu's thesis Section 4.8
% Here, we model both materials with Stiffened gas equation of state with pstiff = 0.

% units: g, mm, s, Pa

% Fluid materials
% materials = [
%   VarFcnJWL(0.28, 3.712e11, 3.23e9, 4.15, 0.95, 1.63e-3);
%   VarFcnSG(7.15, 3.309e8);
% ];

materials = [
  VarFcnSG(1.4, 0);
  VarFcnSG(1.667, 0);
];

% Initial Conditions
% User inputs
rhoIC = [1.0 1.0];
uIC = [0 0];
pIC = [500 0.2];

% Spatial domain
N = 1000;
x_in = 0; x_out = 1;
dx = 1/N;
x = (x_in + dx/2):dx:(x_out-dx/2);

% Temporal parameters
t      = 0; 
dt     = 1e-5; 
tf     = 0.015;
n      = 0;
nsteps = round(tf/dt);

% Interface location
% NOTE: `x_int` stores the x co-ordinate of the material interface. This example has only one.
% However, more material interfaces can be present. In such cases, `x_int` should contain
% x co-ordinate of all interfaces (i.e. an array).
x_interface = [0.5];
num_of_interfaces = length(x_interface);

% Check not needed but added for completeness
if numel(materials) ~= num_of_interfaces+1
    if numel(materials) < num_of_interfaces+1
        error("*** Error: Required fluid materials are %d, whereas only %d were provided.", ...
            num_of_interfaces+1, numel(materials));
    elseif numel(materials) > num_of_interfaces+1
        error("*** Error: Required interface locations are %d, whereas only %d were provided.", ...
            numel(materials)-1, num_of_interfaces);
    end 
end

start_time = tic;

% Setup initial state matrix
% Only 3 primitive variables; [\rho, u_x, p]
Vn = zeros(3, N);
% Not used. Kept here for reference.
% Only 3 conservative variables; [\rho, \rho u_x, E]
% Un = zeros(3, N);
% Stores material ids
Idn = zeros(1, N);

AllIds = 1:numel(materials);

% Populate primitive variables and material ids
Vn(1, :) = rhoIC(end)*ones(1, N);
Vn(2, :) = uIC(end)*ones(1, N);
Vn(3, :) = pIC(end)*ones(1, N);
Idn = AllIds(end)*ones(size(Idn));
index_bounds = ones(1, num_of_interfaces+1);
for j = 1:num_of_interfaces    
    index_bounds(j+1) = find(diff(sign(x - x_interface(j))));
    
    index = index_bounds(j):index_bounds(j+1);
    
    Vn(1, index) = rhoIC(j);
    Vn(2, index) = uIC(j);
    Vn(3, index) = pIC(j);
    Idn(index) = AllIds(j);
end

% Interface x co-ordinate -- might be redundant; could just use Idn
Xn = x_interface;

% For post-processing
% indices for accessing state variables from post processing matrices
keys = reshape(1:3*N, 3, N);
VVn = zeros(3*N, nsteps+1);
XXn = zeros(num_of_interfaces, nsteps+1);

VVn(:, 1) = reshape(Vn, [3*N, 1]);
XXn(:, 1) = reshape(Xn, [num_of_interfaces, 1]);

% =================================================================================================
% Time loop
% =================================================================================================
while t<=tf

    n = n + 1;
    t = t + dt;
    
    fprintf("Step %d: t = %e, dt = %e. Computation time: %.4e s.\n", n, t, dt, toc(start_time));
    
    u = Vn(2, :);
    p = Vn(3, :);
    
%     Interface Tracking Routine 
    Xnp = Xn;
    Idnp = Idn;

%     Advect the material interface(s)
    for j = 1:num_of_interfaces
        
        Xn(j) = real(Xn(j));
        index = find(diff(sign(x - Xn(j))));
        
        x_low  = x(index);
        x_high = x(index+1);
        u_low  = u(index);
        u_high = u(index+1);
        
        c1 = (Xn(j) - x_low)/(x_high - x_low);
        c2 = 1 - c1;
        
        u_interface = c1*u_high + c2*u_low;
        Xnp(j) = Xn(j) + u_interface*dt;
        
    end
    
%     Update material ids
    index_bounds = ones(1, num_of_interfaces+1);
    for j = 1:num_of_interfaces    
        index_bounds(j+1) = find(diff(sign(x - Xnp(j))));

        index = index_bounds(j):index_bounds(j+1);
        Idnp(index) = AllIds(j);
    end
    
%     Phase update -- refer Arthur Rallu's thesis.
    for i=1:N
       if Idnp(i) ~= Idn(i)
%           Compute state from neighbors
           sum_weight    = 0.0;
           sum_weight_Vn = zeros(3, 1);
           
%            This part of the code is taken directly from M2C.
           for neighi = (i-1):(i+1)
               if Idnp(neighi) ~= Idnp(i)
                   continue;
               end
               
               if Idnp(neighi) ~= Idn(neighi)
%                    this node also changed ID skip. This also skips [i]
                   continue;
               end
               
%                Upwinding
               u_neigh = Vn(2, neighi);
               u_neigh_norm = norm(u_neigh);

               if u_neigh_norm ~= 0
                   u_neigh = u_neigh/u_neigh_norm;
               end

               xneighxi = x(i) - x(neighi);
               xneighxi = xneighxi/norm(xneighxi);
               weight   = max(0.0, xneighxi*u_neigh);

               if weight > 0
                   sum_weight = sum_weight + weight;
                   sum_weight_Vn = weight*Vn(:, neighi);
               end
               
           end
           
           if sum_weight == 0
               error("*** Error: Unable to update phase change at (%d)(%e) by " ...
                   + "extrapolation w/ upwinding.", i, x(i));
           end
           Vn(:, i) = sum_weight_Vn/sum_weight;
       end
    end
    
%     Advance in time (uses forward euler method and Exact Riemann solver)
    Vn  = AdvanceOneTimeStep(Vn, materials, x, Idnp, dt);
    Idn = Idnp;
    Xn  = Xnp;
    
%     Save data and update loop
    VVn(:, n) = reshape(Vn, [3*N, 1]);
    XXn(:, n) = reshape(Xn, [num_of_interfaces, 1]);

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
for k = 1:num_of_interfaces
    xline(XXn(k, end),'--', 'linewidth', 1)
end
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
for k = 1:num_of_interfaces
    xline(XXn(k, end),'--', 'linewidth', 1)
end
hold off;
title('Pressure')
txt = ['t = ' num2str(t_vec(end))];
text(0.1,0.1,txt, 'Units', 'normalized');

% Velocity plot
figure(2)
grid on;
box on;
hold on;
plot(x, ut(:, end), 'LineWidth', 2);
for k = 1:num_of_interfaces
    xline(XXn(k, end),'--', 'linewidth', 1)
end
hold off;
title('Velocity')
txt = ['t = ' num2str(t_vec(end))];
text(0.1,0.1,txt, 'Units', 'normalized');

% Post processing (animation)
close all
ff = 20; % frames per second
for n = 1:ff:nsteps
% % % Toggle comments to plot different variables % % %
%     plot(x, rhot(:, n), 'linewidth', 2); a = 1;
%     plot(x, ut(:, n), 'linewidth', 2); a = 2;
    plot(x, pt(:, n), 'linewidth', 2); a = 3;
    for k = 1:num_of_interfaces
        xline(XXn(k, n),'--', 'linewidth', 1)
    end
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