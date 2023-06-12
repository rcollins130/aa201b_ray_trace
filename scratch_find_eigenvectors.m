%%
% FINAL PROJECT: AIRCRAFT NOISE RAY TRACING 3D
% AA201B STANFORD SPRING 2023
% ROBERT COLLINS

% setup
clear;
%close all;

% sim parameters
dt = 0.1; %s
t_max = 5; %s
dz = 0.1; %m
z_range = [0, 1000];

% source and receiver
x_src = [0     0   300];
x_rcv = [50,100,1];
rad_rcv = 15;

% wind
wind_azi = deg2rad(0);

% setup environment
nt = t_max / dt + 1;
z = linspace(z_range(1), z_range(2), (diff(z_range))/dz + 1)';

[T,c,v] = generate_profiles(z, wind_azi);
v = v*0;

% initial ray sweep
n_rays_ele = 7;
n_rays_azi = 8;

min_ele = -pi/2;
max_ele = pi/2;
min_azi = 0;
max_azi = 2*pi - 2*pi/n_rays_azi;

best_d = 9e9;

figure(2); clf; hold on

while best_d > rad_rcv

la_ele = linspace(min_ele, max_ele, n_rays_ele);
la_azi = linspace(min_azi, max_azi, n_rays_azi);

la = zeros(n_rays_ele*n_rays_azi, 2);
ele_idx = zeros(n_rays_ele*n_rays_azi, 1);
azi_idx = zeros(n_rays_ele*n_rays_azi, 1);

ii_ray = 1;
for ii_theta = 1:n_rays_ele
    for ii_phi = 1:n_rays_azi
        la(ii_ray,:) = [la_ele(ii_theta), la_azi(ii_phi)];
        ele_idx(ii_ray) = ii_theta;
        azi_idx(ii_ray) = ii_phi;
        ii_ray = ii_ray +1;
    end
end

% sweep rays
r = nan(nt, 3, size(la, 1));
max_ii_t = repmat(nt, size(la,1),1);

for ii_ray = 1:size(la,1)
    [r(:,:,ii_ray), max_ii_t(ii_ray)]= raytrace_rac( ...
        dt, nt, dz, x_src, la(ii_ray,1), la(ii_ray,2), z, v, c ...
        );
end

% plot initial rays
% Plot Rays
cmap = parula(n_rays_ele);
ax = gca(); 
axis equal
%view([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')
colormap(cmap)
cb = colorbar();
clim([-pi/2, pi/2])
cb.Ticks = -pi/2:pi/8:pi/2;
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
plot3(x_src(1), x_src(2), x_src(3),'ro');
plot3(x_rcv(1), x_rcv(2), x_rcv(3),'bo');

% find ray that got closest to the point
track_dists = vecnorm(x_rcv - r, 2,2);
d_mins = min(track_dists,[],1);
[best_d, best_ray] = min(d_mins)

plot3(ax, ...
r(1:max_ii_t(best_ray), 1, best_ray),...
r(1:max_ii_t(best_ray), 2, best_ray),...
abs(r(1:max_ii_t(best_ray), 3, best_ray)),...
'Color', 'red'...
);

% get neighbors of best ray
next_ele_idx = mod(ele_idx(best_ray) + [-1,1] -1, length(la_ele))+1;
next_azi_idx = mod(azi_idx(best_ray) + [-1,1] -1, length(la_azi))+1;
next_ele = la_ele(next_ele_idx);
next_azi = la_azi(next_azi_idx);
if next_azi(1) > next_azi(2)
    next_azi(1) = next_azi(1) - 2*pi;
end

min_azi= next_azi(1);
max_azi = next_azi(2);
min_ele = next_ele(1);
max_ele = next_ele(2);

n_rays_ele = 5;
n_rays_azi = 5;
end

eig_v = r(1:max_ii_t(best_ray), 1, best_ray);
