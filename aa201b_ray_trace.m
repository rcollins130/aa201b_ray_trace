%%
% FINAL PROJECT: AIRCRAFT NOISE RAY TRACING
% AA201B STANFORD SPRING 2023
% ROBERT COLLINS

%% SETUP
% clear; close all;

%% PHYSICAL PARAMETERS

% Physical Constants
%   Schafer and Vorlander
% temperature gradient
T_0 = 288.15; % Reference sea level temperature, K
L_b = 0.0065; % Temperature lapse rate, K/m

gamma= 1.4; % adiabatic index
R = 287.058; % molar gas constant, dry air, m2/s2

% velocity gradient
% https://en.wikipedia.org/wiki/Log_wind_profile
% https://en.wikipedia.org/wiki/Wind_profile_power_law

% Karman Constant
kar = 0.40; 

% Surface roughness, m 
% https://en.wikipedia.org/wiki/Roughness_length
z_0 = 0.1; 

% Friction velocity
% made this up for now
v_s = 0.05;

%% SIMULATION PARAMETERS
% x, z, t steps
dx = 100; %m
dz = 100; %m
dt = 0.1; %s

% x, z max values
% note x is symmetric about 0
z_max = 10000; %m 
x_max = 5000; %m
t_max = 1000; %s

% source location relative to receiver
source_loc = [-500, 5000]; %[x,z], m

% 


%% DEFINE ENVIRONMENT
% x and z spaces
z = 0:dz:z_max;
% TODO: do I need an X space?
x = -x_max:dx:x_max; % TODO: might cause issues with odd dx/x_max vals

t = 0:dt:t_max;

[X,Z] = meshgrid(x,z);

% temperature profile
T = T_0 - L_b * z;

% velocity profile
% TODO: ADD POWER LAW ABOVE ~1000M?
V = v_s/kar .* log(z/z_0);
V(z<z_0) = 0;

% pressure profile
%   https://en.wikipedia.org/wiki/Barometric_formula
% TODO: check this works with the temp. profile defined above
z_b = 0; % reference height, m
a = 2.2558; % exponent
L_b = 0.0065; % temperature lapse rate, K/m
T_b = 288.15; % sea level reference temperature, K
P_b = 101325; % static pressure, Pa
P = P_b .* ((T_b - (z-z_b).*L_b)).^a;

% speed of sound profile
C = sqrt(gamma*R * T);

% source location
[~,src_row] = min(abs(z-source_loc(2)));
[~,src_col] = min(abs(x-source_loc(1)));

% receiver location
rcv_row = 0;
rcv_col = x_max/dx + 1;

% source level

% define sim parameters

% define plotting parameters

%% INITIAL PLOTTING
% Plot environmental profile
figure(1)
subplot(1,4,1)
plot(T,z)
xlabel('Temperature, K')
ylabel('Altitude, m')

subplot(1,4,2)
plot(P,z)
xlabel('Pressure, Pa')
ylabel('Altitude, m')

subplot(1,4,3)
plot(C,z)
xlabel('Sound Speed, m/s')
ylabel('Altitude, m')

subplot(1,4,4)
plot(V,z)
xlabel('Wind Velocity, m/s')
ylabel('Altitude, m')


%% RAY TRACE
% to start, maybe just iterate thru z vals first?
% initial ray vector directions
n_rays = 5;
theta = linspace(5*pi/4, 7*pi/4, n_rays);

% initialze ray locations & vectors
ray_r = zeros(length(t),2,n_rays);
ray_dr = zeros(size(ray_r));

% initial values
ray_r(1,:,:) = repmat(source_loc,1,1,n_rays);
ray_dr(1,:,:) = reshape([-dz ./ tan(theta); repmat(-dz, 1, length(theta))],[1,2,length(theta)]);

% iterate thru time
for ii_t = 2:length(t)
    % propagate using last dr
    ray_r(ii_t,:,:) = ray_r(ii_t-1,:,:) + ray_dr(ii_t-1,:,:);

    % update dr based on current wavefront position
    ray_dr(ii_t,:,:) = ray_dr(ii_t-1,:,:)-sqrt(ii_t)/1000;
end

figure(2); clf; hold on
for ii_ray = 1:n_rays
    plot(ray_r(:,1,ii_ray), ray_r(:,2,ii_ray))
end
