%%
% FINAL PROJECT: AIRCRAFT NOISE RAY TRACING 3D
% AA201B STANFORD SPRING 2023
% ROBERT COLLINS

%% SETUP
clear; 
% close all;

%% SIMULATION PARAMETERS
dz = 1; % m
dt = 0.01; %s

z_range = [0, 3500]; %m
x_range = [-5000, 5000]; %m
y_range = [-5000, 5000];
t_max = 600; %s

% 
n_rays_ele = 50;
n_rays_azi = 1;

min_ele_deg = -90;
max_ele_deg = 90;

min_ele = deg2rad(min_ele_deg);
max_ele = deg2rad(max_ele_deg);
min_azi = 0;
max_azi = 0;

use_rk4 = 1;

% source location
r_0 = [0,0,500];

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

% Friction velocity, m/s
%   Schafer and Vorlander
v_s = 0.6;

%% DEFINE ENVIRONMENT
z = linspace(z_range(1), z_range(2), (diff(z_range))/dz + 1)';
t = linspace(0, t_max, (t_max)/dt + 1)';

round2zlim = @(ii_z)(min(max(ii_z,1), size(z,1)));

% temperature profile
T = T_0 - L_b * z;

% velocity profile
% TODO: ADD POWER LAW ABOVE ~1000M?
vx= v_s/kar .* log(z/z_0);
vx(z<z_0) = 0;

v = [vx, zeros(size(vx))]; % x, y

% sound speed profile
c = sqrt(gamma * R * T);
% c = -sawtooth(2*pi*z/max(z), 1/2)*10 + 350;
% c = 343*ones(size(z));

% compute derivatives of profiles
v_dz = [diff(v,1,1)/dz; zeros(1,2)];
c_dz = [diff(c)/dz; 0];

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
% ray angle rel2 positive x
la_ele = linspace(min_ele, max_ele, n_rays_ele);
la_azi = linspace(min_azi, max_azi, n_rays_azi);

n_rays = n_rays_ele * n_rays_azi;
la = zeros(n_rays_ele*n_rays_azi, 4);
ii_ray = 1;
for ii_theta = 1:n_rays_ele
    for ii_phi = 1:n_rays_azi
        la(ii_ray,:) = [la_ele(ii_theta), la_azi(ii_phi), ii_theta, ii_phi];
        ii_ray = ii_ray +1;
    end
end

r = zeros(size(t,1), 3, n_rays);
max_ii_t = repmat(size(t,1), n_rays,1);
%r(1,:,:,:) = repmat(r_0, 1, 1, n_rays_theta, n_rays_phi);

tic 
for ii_ray = 1:n_rays    
    [r(:,:,ii_ray), max_ii_t(ii_ray)]= raytrace_rac( ...
        dt, size(t,1), dz, r_0, la(ii_ray,1), la(ii_ray,2), z, v, c ...
        );
end
toc

cmap = parula(n_rays_ele);
figure(2); clf; hold on
ax = gca();
axis equal
view([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')
colormap(cmap)
cb = colorbar();
dele = (max_ele - min_ele)/(n_rays_ele-1);
clim([min_ele, max_ele])
% cb.Ticks = la_ele;
cb.TickLabels = compose("%0.2fº", rad2deg(cb.Ticks));
ylabel(cb, 'Ray Launch Elevation')

for ii_ray = 1:n_rays
    plot3(ax, ...
            r(1:max_ii_t(ii_ray), 1, ii_ray),...
            r(1:max_ii_t(ii_ray), 2, ii_ray),...
            abs(r(1:max_ii_t(ii_ray), 3, ii_ray)),...
            'Color', cmap(la(ii_ray,3),:)...
            );
        pause(0.01)
end
view([0,-1,0])
