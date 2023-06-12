%%
% Validate ray tracer against other model

%% SIMULATION PARAMETERS
% time steps
dt = 0.01; %s
t_max = 10; %s
nt = t_max / dt + 1;

% profiles
dz = 0.001; % m
z_range = [0, 3500]; %m

% rays
n_rays_ele = 20;
n_rays_azi = 1;
min_ele = deg2rad(-90);
max_ele = deg2rad(90);
min_azi = deg2rad(0);
max_azi = deg2rad(0);

% source location
r_0 = [0,0,200];

% wind azimuth
wind_azi_sweep = deg2rad([0, 180, 90]);

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

[T,c] = generate_profiles(z);

% setup launch angles
[la, ele_idx, ~] = launch_angles(...
    min_ele, max_ele, n_rays_ele, ...
    min_azi, max_azi, n_rays_azi);

figure(1); clf;
cmap = parula(n_rays_ele);

for ii_sweep = 1:length(wind_azi_sweep)
    wind_azi = wind_azi_sweep(ii_sweep);
    % velocity profile
    % TODO: ADD POWER LAW ABOVE ~1000M?
    v_mag= v_s/kar .* log(z/z_0);
    v = v_mag * [cos(wind_azi), sin(wind_azi)]; % x, y
    v(z<z_0, :) = 0;
    
    %% MY RAY TRACE 
    r_mine = zeros(nt, 3, size(la,1));
    max_ii_t = repmat(nt, size(la,1),1);
    
    tic 
    for ii_ray = 1:size(la,1)    
        [r_mine(:,:,ii_ray), max_ii_t(ii_ray)]= raytrace_rac( ...
            dt, nt, dz, r_0, la(ii_ray,1), la(ii_ray,2), z, v, c ...
            );
    end
    toc
    
    %% THEIR RAY TRACE 
    r_wilson = zeros(nt, 3, size(la,1));
    tic
    for ii_ray = 1:size(la,1)
        [r_wilson(:,1,ii_ray), ...
        r_wilson(:,2,ii_ray), ...
        r_wilson(:,3,ii_ray)] = raytrace(...
            0:dt:t_max, la(ii_ray,2), la(ii_ray,1), r_0(3), z', c', v(:,1)', v(:,2)',0 ...
            );
    end
    toc
    
    %% SHOWDOWN
    subplot(2,2,ii_sweep); hold on;
    for ii_ray = 1:size(la,1)
        plot(r_mine(1:max_ii_t(ii_ray),1,ii_ray), abs(r_mine(1:max_ii_t(ii_ray),3,ii_ray)), ...
            '-', 'LineWidth',2, 'Color', cmap(ele_idx(ii_ray),:));
        plot(r_wilson(:,1,ii_ray), abs(r_wilson(:,3,ii_ray)), 'r--');
    end
    xlabel('x (m)')
    ylabel('z (m)')
    colormap(cmap)
    cb = colorbar;
    clim([min_ele, max_ele])
    cb.Ticks = min_ele:pi/4:max_ele;
    cb.TickLabels = compose("%0.1f",rad2deg(cb.Ticks));
    ylabel(cb, 'Launch Elevation')
    axis equal
    grid on
    ylim(z_range)
    xlim([0,3500])
end

subplot(2,2,4); hold on;
for ii_ray = 1:size(la,1)
    plot(r_mine(1:max_ii_t(ii_ray),1,ii_ray), abs(r_mine(1:max_ii_t(ii_ray),2,ii_ray)), ...
        '-', 'LineWidth',2, 'Color', cmap(ele_idx(ii_ray),:));
    plot(r_wilson(:,1,ii_ray), abs(r_wilson(:,2,ii_ray)), 'r--');
end
xlabel('x (m)')
ylabel('y (m)')
colormap(cmap)
cb = colorbar;
clim([min_ele, max_ele])
cb.Ticks = min_ele:pi/4:max_ele;
cb.TickLabels = compose("%0.1f",rad2deg(cb.Ticks));
ylabel(cb, 'Launch Elevation')
grid on
xlim([0,3500])
