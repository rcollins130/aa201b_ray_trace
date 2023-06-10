%%
% FINAL PROJECT: AIRCRAFT NOISE RAY TRACING SNELL IMPLEMENTATION
% AA201B STANFORD SPRING 2023
% ROBERT COLLINS

%% SETUP
clear; %close all;

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
dz = 10; %m

% x, z max values
% note x is symmetric about 0
z_max = 10000; %m 
x_max = 5000; %m
t_max = 1000; %s

% source location relative to receiver
source_loc = [0, 7500]; %[x,z], m

% 

%% DEFINE ENVIRONMENT
% define z space 
z = (0:dz:z_max)';

% map source location to z index
[~, src_ii_z] = min(abs(source_loc(2)-z));

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
% C = sqrt(gamma*R * T);
%C = ones(size(z)) * 350;
C = -sawtooth(2*pi*z/max(z), 1/2)*50 + 350;
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
n_rays = 100;
max_trace = 1e3;
theta = linspace(0, pi, n_rays);

% initialze ray locations & vectors
% [pidx, (t, x, z, dx, dz, ii_z, theta), ray]
ray_data = zeros(max_trace,6,n_rays);

% initial values
ray_data(1,:,:) = [
    repmat([0, source_loc],1,1,n_rays) ...
    reshape(-sign(abs(theta')-pi/2) .* dz .* tan(theta'), [1,1,n_rays]) ...
    reshape(sign(abs(theta')-pi/2) .* dz, [1,1, n_rays]) ...
    repmat(src_ii_z, 1, 1, n_rays)
    ];

figure(2); clf; hold on;
for ii_ray = 1:n_rays
    ray_lines(ii_ray) = plot(ray_data(:,2,ii_ray), ray_data(:,3,ii_ray));
end
axis equal
ylim([0,z_max])
xlim([-50000, 50000])

% iterate thru each ray
for ii_ray = 1:n_rays
    % trace ray
    for ii_pt = 2:max_trace
        % extract last ray
        last_ray = ray_data(ii_pt-1, :, ii_ray);

        % calculate new position  
        this_x = last_ray(2) + last_ray(4);
        this_ii_z = last_ray(6) + sign(last_ray(5)); 
        this_z = z(this_ii_z);

        % calculate new t
        dt = vecnorm([this_x, this_z] - last_ray(2:3)) / C(this_ii_z);
        this_t = last_ray(1) + dt;
        
        % apply snell's law
        this_theta = asin(-sign(last_ray(5)) * last_ray(4) ./ norm(last_ray(4:5)));
        next_sin_theta = C(this_ii_z)/C(this_ii_z-sign(last_ray(5))) * sin(this_theta);
        if abs(next_sin_theta) >1
            this_dz = -last_ray(5);
            this_dx = this_dz * tan(this_theta);
        else
            this_dz = last_ray(5);
            this_dx = -this_dz * tan(asin(next_sin_theta));
        end

        % update dr based on current wavefront position
        

        % update next data 
        ray_data(ii_pt, :, ii_ray) = [
            this_t, this_x, this_z, this_dx, this_dz, this_ii_z
            ];
        if this_ii_z <= 2 || this_ii_z >= size(z,1)-1
            break
        end
        
        % some debugging help
        tracehelp = ray_data(1:ii_pt,:,ii_ray);
        ray_lines(ii_ray).XData = real(ray_data(:,2,ii_ray));
        ray_lines(ii_ray).YData = real(ray_data(:,3,ii_ray));
        pause(0.001);
    end
end

Xs = reshape(ray_data(:,2,:), [size(ray_data,1)*size(ray_data,3),1]);
Zs = reshape(ray_data(:,3,:), [size(ray_data,1)*size(ray_data,3),1]);
Ts = reshape(ray_data(:,1,:), [size(ray_data,1)*size(ray_data,3),1]);

[Xq,Zq] = meshgrid(linspace(-50000, 50000,100), z);
T = griddata(Xs,Zs,Ts,Xq,Zq);
figure
contourf(Xq,Zq,T)
axis equal

% figure(2); clf; hold on
% for ii_ray = 1:n_rays
%     plot(ray_data(:,2,ii_ray), ray_data(:,3,ii_ray))
% end
% 
% ylim([0,z_max])
% xlim([-15000, 15000])