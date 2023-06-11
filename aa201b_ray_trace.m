%%
% FINAL PROJECT: AIRCRAFT NOISE RAY TRACING
% AA201B STANFORD SPRING 2023
% ROBERT COLLINS

%% SETUP
clear; 
% close all;
warning('this code does not work!')
%% SIMULATION PARAMETERS
dz = 10; % m
dt = 0.01; %s

z_range = [0, 10000]; %m
x_range = [-1000, 1000]; %m
t_max = 60; %s

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

%% DEFINE ENVIRONMENT
z = linspace(z_range(1), z_range(2), (diff(z_range))/dz + 1)';
t = linspace(0, t_max, (t_max)/dt + 1)';

% temperature profile
T = T_0 - L_b * z;

% velocity profile
% TODO: ADD POWER LAW ABOVE ~1000M?
vx= v_s/kar .* log(z/z_0)* 0;
vx(z<z_0) = 0;

% sound speed profile
% c = sqrt(gamma * R * T);
c = -sawtooth(2*pi*z/max(z), 1/2)*50 + 350;

% compute derivatives of profiles
vx_dz = [0; diff(vx)];
c_dz = [0; diff(c)];

%% INITIAL PLOTTING
% Plot environmental profile
figure(1); clf;
subplot(1,4,1)
plot(T,z)
xlabel('Temperature, K')
ylabel('Altitude, m')

subplot(1,4,3)
plot(c,z)
xlabel('Sound Speed, m/s')
ylabel('Altitude, m')

subplot(1,4,4)
plot(vx,z)
xlabel('Horizontal Wind Velocity, m/s')
ylabel('Altitude, m')

%% RAY TRACE

% 
n_rays = 150;

% ray angle rel2 positive x
launch_angle = linspace(pi, 2*pi, n_rays);

figure(2); clf; hold on

for ii_ray = 1:n_rays

% initialize ray and normal
r = zeros(size(t, 1), 2);
n = zeros(size(t, 1), 2);

% initial ray position 
r(1,:) = [0,6000]; % x, z
n(1,:) = [cos(launch_angle(ii_ray)), sin(launch_angle(ii_ray))]; % u, w

for ii_t = 1:size(t,1)-1
    % get ii_z 
    [~, ii_z] = min(abs(z-r(ii_t, 2)));
    
    if ii_z == 1 || ii_z == size(z, 1)
        r = r(1:ii_t, :);
        break
    end

    % current values
    s = n(ii_t,:) ./ (c(ii_z) + n(ii_t,:) * vx(ii_z,:)');
    omega = 1 - vx(ii_z) * s(1);
    sz_dt = - omega./c(ii_z) * c_dz(ii_z) - s(1) .* vx_dz(ii_z);
    
    % advance values (euler method)
    r(ii_t+1, :) = r(ii_t,:) + c(ii_z) .* n(ii_t,:) * dt;
    s(2) = s(2) + sz_dt*dt;
    n(ii_t+1, :) = s/norm(s);

    % advance values (RK4)

end

plot(r(:,1), r(:,2))

end