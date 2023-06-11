%%
% FINAL PROJECT: AIRCRAFT NOISE RAY TRACING 3D
% AA201B STANFORD SPRING 2023
% ROBERT COLLINS

%% SETUP
clear; 
% close all;

%% SIMULATION PARAMETERS
dz = 0.1; % m
dt = 0.1; %s

z_range = [0, 10000]; %m
x_range = [-1000, 1000]; %m
y_range = [-1000, 1000];
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
vx= v_s/kar .* log(z/z_0)* 10;
vx(z<z_0) = 0;

v = [vx, zeros(size(vx)), zeros(size(vx))]; % x, y, z

% sound speed profile
%c = sqrt(gamma * R * T);
%c = -sawtooth(2*pi*z/max(z), 1/2)*10 + 350;
c = 343*ones(size(z));

% compute derivatives of profiles
v_dz = [zeros(1,3); diff(v,1,1)/dz];
c_dz = [0; diff(c)/dz];

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
n_rays_theta = 15;
n_rays_phi = 15;

min_theta = 17/32*pi;
max_theta = 3/4*pi;
min_phi = 0;
max_phi = 2*pi;

% ray angle rel2 positive x
la_theta = linspace(min_theta, max_theta, n_rays_theta);
la_phi = linspace(min_phi, max_phi, n_rays_phi+1);
la_phi = la_phi(1:end-1);

r_0 = [0,0,5000];

cmap = parula(n_rays_theta);
figure(2); clf; hold on
ax = gca();
axis equal
view([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')
colormap(cmap)
cb = colorbar();
dtheta = (max_theta - min_theta)/(n_rays_theta-1);
clim([min_theta-dtheta/2, max_theta+dtheta/2])
cb.Ticks = la_theta;
cb.TickLabels = compose("%0.2f \\pi", la_theta/pi);
ylabel(cb, 'Ray Launch Elevation')
r = zeros(size(t,1), 3, n_rays_theta, n_rays_phi);
r(1,:,:,:) = repmat(r_0, 1, 1, n_rays_theta, n_rays_phi);

% plot_array = cell(n_rays_theta, n_rays_phi);
% for ii_phi = 1:n_rays_phi
% for ii_theta = 1:n_rays_theta
%     plot_array{ii_theta, ii_phi} = plot3( ...
%         r(:,1,ii_theta,ii_phi), ...
%         r(:,2,ii_theta,ii_phi), ...
%         r(:,3,ii_theta,ii_phi), ...
%         'Color', cmap(ii_theta, :) ...
%     );
% end
% end

for ii_phi = 1:n_rays_phi
for ii_theta = 1:n_rays_theta

    % initialize ray vector
    % r = zeros(size(t, 1), 3);
    % n = zeros(size(t, 1), 3);

    % initial ray normal 
    n = [
        sin(la_theta(ii_theta))*cos(la_phi(ii_phi)) ...
        sin(la_theta(ii_theta))*sin(la_phi(ii_phi)) ...
        cos(la_theta(ii_theta))
    ]; % x, y ,z

    % compute starting ii_z
    [~, ii_z] = min(abs(z-r(1, 3)));

    % compute initial slowness vector
    s = n ./ (c(ii_z) + n * v(ii_z,:)');
    
    % compute slowness perp and omega vec
    s_perp = s(1:2);
    s_z = s(3);
    omega = 1 - v(:,1:2) * s(1:2)';
    
    % compute slowness z derivative
    sz_dt = - omega ./ c .* c_dz - v_dz(:,1:2) * s_perp'; 
    
    max_ii_t = size(t,1);

    % iterate thru time
    for ii_t = 1:size(t,1)-1
        % get current ii_z 
        [~, ii_z] = min(abs(z-r(ii_t, 3, ii_theta, ii_phi)));
        
        if ii_z == 1 || ii_z == size(z, 1)
            max_ii_t = ii_t;
            break
        end

        % EULER METHOD
        % current values
        % sz_dt = - omega(ii_z)./c(ii_z) * c_dz(ii_z) - s(1:2) * v_dz(ii_z,1:2)';
        
        % advance values (euler method)
        r(ii_t+1, :, ii_theta, ii_phi) = r(ii_t,:, ii_theta, ii_phi) + ...
            (c(ii_z).^2 / omega(ii_z) .* [s_perp, s_z] + v(ii_z,:)) .* dt;
        s_z = s_z + sz_dt(ii_z) * dt;

        % advance values (RK4)
        
        % compute next s_z
        % z_1 = z(ii_z);
        % a_1 = g(ii_z);
        % 
        % z_2 = z_1 + dt* a_1/2;
        % [~, ii_z2] = min(abs(z-z_2));
        % a_2 = g(ii_z2);
        % 
        % z_3 = z_1 + dt* a_2/2;
        % [~, ii_z3] = min(abs(z-z_3));
        % a_3 = g(ii_z3);   
        % 
        % z_4 = z_1 + dt* a_3;
        % [~, ii_z4] = min(abs(z-z_4));
        % a_4 = g(ii_z4);
        % 
        % s_next = s;
        % s_next(3) = s(3) + dt/6 * (a_1 + 2*a_2 + 2*a_3 + a_4);
        % s_h = 1/2*(s+s_next);
        % 
        % % compute next r
        % b_1 = c(ii_z).^2 / omega(ii_z) * s + v(ii_z,:);
        % 
        % z_2 = z_1 + dt* b_1(3)/2;
        % [~, ii_z2] = min(abs(z-z_2));
        % b_2 = c(ii_z2).^2 / omega(ii_z2) * s_h + v(ii_z2,:);
        % 
        % z_2 = z_1 + dt* b_2(3)/2;
        % [~, ii_z3] = min(abs(z-z_2));
        % b_3 = c(ii_z3).^2 / omega(ii_z3) * s_h + v(ii_z3,:);
        % 
        % z_4 = z_1 + dt* b_3(3);
        % [~, ii_z4] = min(abs(z-z_4));
        % b_4 = c(ii_z4).^2 / omega(ii_z4) * s_next + v(ii_z4,:);
        % 
        % r_next = r(ii_t,:) + dt/6 * (b_1 + 2*b_2 + 2*b_3 + b_4);

        % % advance values
        % s = s_next;
        % r(ii_t+1, :) = r_next;

        % update plot
        % plot_array{ii_theta, ii_phi}.XData(ii_t) = r(ii_t, 1, ii_theta, ii_phi);
        % plot_array{ii_theta, ii_phi}.YData(ii_t) = r(ii_t, 2, ii_theta, ii_phi);
        % plot_array{ii_theta, ii_phi}.ZData(ii_t) = r(ii_t, 3, ii_theta, ii_phi);
        % pause(0.0001)
    end
    
    plot3( ...
        r(1:max_ii_t, 1, ii_theta, ii_phi),...
        r(1:max_ii_t, 2, ii_theta, ii_phi),...
        r(1:max_ii_t, 3, ii_theta, ii_phi),...
        'Color', cmap(ii_theta,:)...
        );
    pause(0.01)
end
end

