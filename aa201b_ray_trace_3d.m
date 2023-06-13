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
n_rays_azi = 25;
min_ele = deg2rad(-45);
max_ele = deg2rad(10);
min_azi = deg2rad(0);
max_azi = deg2rad(360-17/360);

% source location
r_0 = [0,0,750];

% wind azimuth
wind_azi = deg2rad(0);

%% DEFINE ENVIRONMENT
% height profile
z = linspace(z_range(1), z_range(2), (diff(z_range))/dz + 1)';

% default environment profiles
[T,c,v] = generate_profiles(z, wind_azi);

% alternate profiles
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
subplot(1,3,1)
plot(T,z)
xlabel('Temperature, K')
ylabel('Altitude, m')

subplot(1,3,2)
plot(c,z)
xlabel('Sound Speed, m/s')
ylabel('Altitude, m')

subplot(1,3,3)
plot(vecnorm(v,2,2),z)
xlabel('Wind Velocity, m/s')
ylabel('Altitude, m')

%% Plot Rays
cmap = parula(n_rays_ele);
figure(2); clf; hold on
%subplot(2,1,2)
ax = gca(); cla; hold on;
axis equal
view([0,-1,0])
xlim([-5000,5000])
ylim([-5000,5000])
zlim([0,3500])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
colormap(cmap)
cb = colorbar();
dele = (max_ele - min_ele)/(n_rays_ele-1);
clim([min_ele, max_ele])
cb.Ticks = -pi:pi/8:pi;
cb.TickLabels = compose("%0.2fº", rad2deg(cb.Ticks));
ylabel(cb, 'Ray Launch Elevation')

for ii_ray = 1:size(la,1)
    msk = r(:,3,ii_ray)>0;
    plot3(ax, ...
            r(msk, 1, ii_ray),...
            r(msk, 2, ii_ray),...
            abs(r(msk, 3, ii_ray)),...
            'Color', cmap(ele_idx(ii_ray),:)...
            );
        pause(0.01)
end
