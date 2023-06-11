%%
% FINAL PROJECT: AIRCRAFT NOISE RAY TRACING 3D
% AA201B STANFORD SPRING 2023
% ROBERT COLLINS

%% SETUP
clear; 
% close all;

%% SIMULATION PARAMETERS
% time steps
dt = 0.01; %s
t_max = 600; %s
nt = t_max / dt + 1;

% profiles
dz = 0.1; % m
z_range = [0, 3500]; %m

% rays
n_rays_ele = 25;
n_rays_azi = 4;
min_ele = deg2rad(-90);
max_ele = deg2rad(90);
min_azi = deg2rad(-180);
max_azi = deg2rad(90);

% source location
r_0 = [0,0,500];

% wind azimuth
wind_azi = deg2rad(0);

%% PHYSICAL PARAMETERS
T_0 = 288.15; % Reference sea level temperature, K
L_b = 0.0065; % Temperature lapse rate, K/m
gamma= 1.4; % adiabatic index
R = 287.058; % molar gas constant, dry air, m2/s2
kar = 0.40; % Karman Constant
z_0 = 0.1; % Surface roughness, m 
v_s = 0.6; % Friction velocity, m/s

%% DEFINE ENVIRONMENT
% height profile
z = linspace(z_range(1), z_range(2), (diff(z_range))/dz + 1)';

% temperature profile
T = T_0 - L_b * z;

% velocity profile
% TODO: ADD POWER LAW ABOVE ~1000M?
v_mag= v_s/kar .* log(z/z_0);
v = -v_mag * [cos(wind_azi), sin(wind_azi)]; % x, y
v(z<z_0, :) = 0;

% sound speed profile
c = sqrt(gamma * R * T);
% c = -sawtooth(2*pi*z/max(z), 1/2)*10 + 350;
% c = 343*ones(size(z));


%% RAY TRACE
[la, ele_idx, ~] = launch_angles(...
    min_ele, max_ele, n_rays_ele, ...
    min_azi, max_azi, n_rays_azi);

r = zeros(nt, 3, size(la,1));
max_ii_t = repmat(nt, size(la,1),1);

tic 
for ii_ray = 1:size(la,1)    
    [r(:,:,ii_ray), max_ii_t(ii_ray)]= raytrace_rac( ...
        dt, nt, dz, r_0, la(ii_ray,1), la(ii_ray,2), z, v, c ...
        );
end
toc

%% PLOTTING
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
plot(vecnorm(v,2,2),z)
xlabel('Wind Velocity, m/s')
ylabel('Altitude, m')

% Plot Rays
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

for ii_ray = 1:size(la,1)
    plot3(ax, ...
            r(1:max_ii_t(ii_ray), 1, ii_ray),...
            r(1:max_ii_t(ii_ray), 2, ii_ray),...
            abs(r(1:max_ii_t(ii_ray), 3, ii_ray)),...
            'Color', cmap(ele_idx(ii_ray),:)...
            );
        pause(0.01)
end
view([0,-1,0])
