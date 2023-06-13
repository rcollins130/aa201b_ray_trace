

%%
clear; close all;

figure(2); clf; hold on

subplot(1,2,2)
ii_wad = 1;
for wad =[0,1]
    subplot(1,2,ii_wad)
    ii_wad = 1+ii_wad;
%% 
dt = 0.01;
t_max = 30; 
nt = t_max / dt+1;
dz = 1;
z_range=[0,5000];

wind_azi = deg2rad(0);
z = linspace(z_range(1), z_range(2), (diff(z_range))/dz + 1)';
[T,c,v] = generate_profiles(z, wind_azi);
v = v*wad;

% observer point
x_src = [0,0,1000];
%x_rcv = [50,2000,2];

% define rays to trace
n_rays_ele = 9;
n_rays_azi = 16;
min_ele = deg2rad(-30);
max_ele = deg2rad(-1);
min_azi = deg2rad(0);
max_azi = deg2rad(360-360/16);

[la, ele_idx, azi_idx] = launch_angles( ...
    min_ele,max_ele, n_rays_ele, ...
    min_azi,max_azi,n_rays_azi);

d_ele = 0.01;
d_azi = 0.01;

%cmap = parula(256);
ax = gca(); cla; hold on;
axis equal
view([0,0,1])
%view([1,1,1])
% xlim([-5000,5000])
% ylim([-5000,5000])
zlim([0,3500])
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
% colormap(cmap)
% cb = colorbar();
% clim([-pi/2, pi/2])
% cb.Ticks = -pi/2:pi/8:pi/2;
% cb.TickLabels = compose("%0.2fº", rad2deg(cb.Ticks));
% ylabel(cb, 'Ray Launch Elevation')

% 
r = nan(nt, 3, size(la,1));
% 
ray_ground_power = nan(size(la,1),3);
for ii_ray=1:size(la,1)
    % trace core ray
    [r(:,:,ii_ray), ~]= raytrace_rac( ...
        dt, nt, dz, x_src, la(ii_ray,1), la(ii_ray,2), z, v, c ...
        );
    
    % find point at ground
    reflected = find(r(:,3,ii_ray)<=0);
    if isempty(reflected)
        continue
    end
    iit_ground = reflected(1);
    xy_ground = r(iit_ground, 1:2, ii_ray);

    plot3(ax, ...
        r(1:iit_ground, 1, ii_ray),...
        r(1:iit_ground, 2, ii_ray),...
        abs(r(1:iit_ground, 3, ii_ray)),...
        'Color', 'k'...%cmap(ceil((la(ii_ray,1)+pi/2)/pi * 255)+1,:)...
        );
    pause(0.01)

    % trace adjacent rays
    r_adj = nan(nt, 3, 3);
    r_adj(:,:,1) = r(:,:,ii_ray);
    ii_adj = 2;
    s_ele = d_ele;
    for s_azi = [-d_azi, d_azi]
        r_adj(:,:,ii_adj) = raytrace_rac( ...
            dt, nt, dz, x_src, la(ii_ray,1)+s_ele, la(ii_ray,2)+s_azi, z, v, c ...
        );  
        
        plot3(ax, ...
            r_adj(1:iit_ground, 1, ii_adj),...
            r_adj(1:iit_ground, 2, ii_adj),...
            abs(r_adj(1:iit_ground, 3, ii_adj)),...
            'k--' ...
        );
        pause(0.01)

        ii_adj = ii_adj+1;
    end
    
    % get area at iit=2;
    xs = squeeze(r_adj(2,1,:))';
    ys = squeeze(r_adj(2,2,:))';
    zs = squeeze(r_adj(2,3,:))';

    ons = [1 1 1];
    A_in = 0.5*sqrt( ...
        det([xs; ys; ons])^2 + ...
        det([ys; zs; ons])^2 + ...
        det([zs; xs; ons])^2 ...
        );

    % get area at iit=t_ground
    xs = squeeze(r_adj(iit_ground,1,:))';
    ys = squeeze(r_adj(iit_ground,2,:))';
    zs = squeeze(r_adj(iit_ground,3,:))';

    ons = [1 1 1];
    A_gr = 0.5*sqrt( ...
        det([xs; ys; ons])^2 + ...
        det([ys; zs; ons])^2 + ...
        det([zs; xs; ons])^2 ...
        );

    ray_ground_power(ii_ray,:) = [xy_ground, A_in / A_gr];
end

rgp = ray_ground_power(~isnan(ray_ground_power(:,3)),:);

[xq,yq] = meshgrid( ...
    linspace(-6000,6000),...
    linspace(-6000,6000));
vq = griddata( ...
    rgp(:,1), ...
    rgp(:,2), ...
    rgp(:,3), ...
    xq,yq);

contourf(xq,yq,10*log(vq))
cb=colorbar;
clim(10*log([1.0026e-07, 2.7851e-06]))
ylabel(cb, "Spreading Loss, dB")
end