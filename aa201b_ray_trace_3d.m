%%
% FINAL PROJECT: AIRCRAFT NOISE RAY TRACING 3D
% AA201B STANFORD SPRING 2023
% ROBERT COLLINS

%% SETUP
clear; 
% close all;

%% SIMULATION PARAMETERS
dz = 0.1; % m
dt = 0.01; %s

z_range = [0, 10000]; %m
x_range = [-1000, 1000]; %m
y_range = [-1000, 1000];
t_max = 500; %s

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
% c = sqrt(gamma * R * T);
c = -sawtooth(2*pi*z/max(z), 1/2)*10 + 350;
%c = 343*ones(size(z));

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
n_rays_theta = 4;
n_rays_phi = 4;

% ray angle rel2 positive x
la_theta = linspace(17/32*pi, 3/4*pi, n_rays_theta);
la_phi = linspace(0, 2*pi, n_rays_phi+1);
la_phi = la_phi(1:end-1);

r_0 = [0,0,5000];

figure(2); clf; hold on
ax = gca();
axis equal
view([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')

for ii_phi = 1:n_rays_phi
for ii_theta = 1:n_rays_theta

    % initialize ray vector
    r = zeros(size(t, 1), 3);
    % n = zeros(size(t, 1), 3);

    % initial ray position 
    r(1,:) = r_0; % x, y, z

    % initial ray normal 
    n = [
        sin(la_theta(ii_theta))*cos(la_phi(ii_phi)) ...
        sin(la_theta(ii_theta))*sin(la_phi(ii_phi)) ...
        cos(la_theta(ii_theta))
    ]; % x, y ,z

    % compute starting ii_z
    [~, ii_z] = min(abs(z-r(1, 3)));
    % 
    % % compute slowness vector
    s = n ./ (c(ii_z) + n * v(ii_z,:)');
    % 
    % % compute slowness perp and omega vec
    % s_perp = s(1:2);
    % omega = 1 - v(:,1:2) * s_perp';
    % 
    % % compute d/dt s_z = g
    % g = - omega ./ c .* c_dz - v_dz(:,1:2) * s_perp'; 

    for ii_t = 1:size(t,1)-1
        % get ii_z 
        [~, ii_z] = min(abs(z-r(ii_t, 3)));
        
        if ii_z == 1 || ii_z == size(z, 1)
            r = r(1:ii_t, :);
            break
        end

        % current values
        % s = n ./ (c(ii_z) + n * v(ii_z,:)');
        omega = 1 - v(ii_z, 1:2) * s(1:2)';
        sz_dt = - omega./c(ii_z) * c_dz(ii_z) - s(1:2) * v_dz(ii_z,1:2)';
        
        % advance values (euler method)
        r(ii_t+1, :) = r(ii_t,:) + (c(ii_z).^2 / omega .* s + v(ii_z,:)) .* dt;
        s(3) = s(3) + sz_dt * dt;
        %n = s/norm(s);

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
    end

    plot3(ax, r(:,1), r(:,2), r(:,3))
    pause(0.01)

end
end

