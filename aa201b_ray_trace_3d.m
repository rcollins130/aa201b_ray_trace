%%
% FINAL PROJECT: AIRCRAFT NOISE RAY TRACING 3D
% AA201B STANFORD SPRING 2023
% ROBERT COLLINS

%% SETUP
clear; 
% close all;

%% SIMULATION PARAMETERS
dz = 1; % m
dt = 0.1; %s

z_range = [0, 3500]; %m
x_range = [-5000, 5000]; %m
y_range = [-5000, 5000];
t_max = 10; %s

% 
n_rays_ele = 21;
n_rays_azi = 1;

min_ele_deg = -85;
max_ele_deg = 85;

min_ele = deg2rad(min_ele_deg);
max_ele = deg2rad(max_ele_deg);
min_azi = 0;
max_azi = 0;

use_rk4 = 1;

% source location
r_0 = [0,0,200];

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

v = [vx, zeros(size(vx)), zeros(size(vx))]; % x, y, z

% sound speed profile
% c = sqrt(gamma * R * T);
%c = -sawtooth(2*pi*z/max(z), 1/2)*10 + 350;
c = 343*ones(size(z));

% compute derivatives of profiles
v_dz = [diff(v,1,1)/dz; zeros(1,3)];
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
    this_r = zeros(size(t, 1), 3);
    this_r(1,:) = r_0;
    % initialize ray vector
    % r = zeros(size(t, 1), 3);
    % n = zeros(size(t, 1), 3);
    
    % initial ray normal 
    % n = [
    %     sin(la(ii_ray, 1)) * cos(la(ii_ray, 2)) ...
    %     sin(la(ii_ray, 1)) * sin(la(ii_ray, 2)) ...
    %     cos(la(ii_ray, 1))
    % ]; % x, y ,z
    
    n = [
        cos(la(ii_ray, 1)) * cos(la(ii_ray, 2)) ...
        -cos(la(ii_ray, 1)) * sin(la(ii_ray, 2)) ...
        sin(la(ii_ray, 1))
    ];

    % compute starting ii_z
    % TODO: possibly move outside
    [~, ii_z0] = min(abs(z-this_r(1, 3)));

    % compute initial slowness vector
    s = n ./ (c(ii_z0) + n * v(ii_z0,:)');
    
    % compute slowness perp and omega vec
    s_perp = s(1:2);
    s_z = s(3);
    omega = 1 - v(:,1:2) * s(1:2)';
    
    % compute slowness z derivative
    sz_dt = - omega ./ c .* c_dz - v_dz(:,1:2) * s_perp'; 
    
    % setup RK4 function
    f = @(sgn, ii_z, s_z) [
        c(ii_z).^2 ./ omega(ii_z) * s_perp(1) + v(ii_z,1); 
        c(ii_z).^2 ./ omega(ii_z) * s_perp(2) + v(ii_z,2); 
        c(ii_z).^2 ./ omega(ii_z) * s_z;
        sgn*(-omega(ii_z)./c(ii_z) * c_dz(ii_z) - s_perp * v_dz(ii_z,1:2)');
    ];
    
    ii_z1 = ii_z0;
    % iterate thru time
    for ii_t = 1:size(t,1)-1
        if ii_z1 >= size(z, 1) || ...
            this_r(ii_t,1) > x_range(2) || ...
            this_r(ii_t,1) < x_range(1) || ...
            this_r(ii_t,2) > y_range(2) || ...
            this_r(ii_t,2) < y_range(1)
            max_ii_t(ii_ray) = ii_t;
            break
        end
        
        if use_rk4
            % advance values (RK4)
            z1 = this_r(ii_t, 3);
            a1 = f(sign(z1), ii_z1, s_z);
    
            z2 = z1 + dt/2 * a1(3);
            ii_z2 = round2zlim(ii_z1 + round((abs(z2)-z1)/dz));
            a2 = f(sign(z2), ii_z2, s_z + dt/2 * a1(4));
    
            z3 = z1 + dt/2 * a2(3);
            ii_z3 = round2zlim(ii_z1 + round((abs(z3)-z1)/dz));
            a3 = f(sign(z3), ii_z3, s_z + dt/2 * a2(4));
    
            z4 = z1 + dt * a3(3);
            ii_z4 = round2zlim(ii_z1 + round((abs(z4)-z1)/dz));
            a4 = f(sign(z4), ii_z4, s_z + dt * a3(4));
    
            dx = dt/6 * (a1 + 2*a2 + 2*a3 + a4);
    
            this_r(ii_t+1, :) = this_r(ii_t, :) + dx(1:3)';
            s_z = s_z + dx(4);
        else
            % advance values (euler method)
            warning('reflections not implemented');
            this_r(ii_t+1, :) = this_r(ii_t,:) + ...
                (c(ii_z1).^2 / omega(ii_z1) .* [s_perp, s_z] + v(ii_z1,:)) .* dt;
            s_z = s_z + sz_dt(ii_z1) * dt;
        end
        
        ii_z1 = round2zlim(ii_z0 + round((abs(this_r(ii_t+1, 3))-r_0(3))/dz));
    end
    r(:,:,ii_ray) = this_r;
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
