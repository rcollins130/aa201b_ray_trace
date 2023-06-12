function [best_r, hit] = find_eigenvector(dt, nt, dz, x_src,x_rcv, rad_rcv, z,v,c)
%FIND_EIGENVECTOR Summary of this function goes here
%   Detailed explanation goes here
max_iters = 15;
do_plot = 0;

% initial sweep pattern
n_rays_ele_init = 7;
n_rays_azi_init = 8;

la_ele = linspace(-pi/2, pi/2, n_rays_ele_init);
la_azi = linspace(0, 2*pi - 2*pi/n_rays_azi_init, n_rays_azi_init);

n_rays_ele = 5;
n_rays_azi = 5;

best_d = 9e9;

ii = 1;
while best_d > rad_rcv && ii<max_iters
    
    % construct rays
    n_rays = length(la_ele)*length(la_azi);
    la = zeros(n_rays, 2);
    ele_idx = zeros(n_rays, 1);
    azi_idx = zeros(n_rays, 1);
    
    ii_ray = 1;
    for ii_ele = 1:length(la_ele)
        for ii_azi = 1:length(la_azi)
            la(ii_ray,:) = [la_ele(ii_ele), la_azi(ii_azi)];
            ele_idx(ii_ray) = ii_ele;
            azi_idx(ii_ray) = ii_azi;
            ii_ray = ii_ray +1;
        end
    end
    
    % trace rays
    r = nan(nt, 3, size(la, 1));
    max_ii_t = repmat(nt, size(la,1),1);
    
    for ii_ray = 1:size(la,1)
        [r(:,:,ii_ray), max_ii_t(ii_ray)]= raytrace_rac( ...
            dt, nt, dz, x_src, la(ii_ray,1), la(ii_ray,2), z, v, c ...
            );
    end

    if do_plot
        cmap = parula(length(la_ele));
        for ii_ray = 1:size(la,1)
        plot3(...
            r(1:max_ii_t(ii_ray), 1, ii_ray),...
            r(1:max_ii_t(ii_ray), 2, ii_ray),...
            abs(r(1:max_ii_t(ii_ray), 3, ii_ray)),...
            'Color', cmap(ele_idx(ii_ray),:)...
            );
        end
    end
    
    % find ray that got closest to the point
    % distances from every point on ray to receiver 
    track_dists = vecnorm(x_rcv - r, 2,2);
    d_mins = min(track_dists,[],1);
    [best_d, best_ray] = min(d_mins);
    
    if do_plot
        plot3(...
        r(1:max_ii_t(best_ray), 1, best_ray),...
        r(1:max_ii_t(best_ray), 2, best_ray),...
        abs(r(1:max_ii_t(best_ray), 3, best_ray)),...
        'Color', 'red'...
        );
    end 

    % get neighbors of best ray
    next_ele_idx = mod(ele_idx(best_ray) + [-1,1] -1, length(la_ele))+1;
    next_azi_idx = mod(azi_idx(best_ray) + [-1,1] -1, length(la_azi))+1;
    next_ele = la_ele(next_ele_idx);
    next_azi = la_azi(next_azi_idx);
    % flip around azimuth if wraps around origin
    if next_azi(1) > next_azi(2)
        next_azi(1) = next_azi(1) - 2*pi;
    end
    % if ele is near nadir, sweep the whole nadir again
    if next_ele(1) > next_ele(2)
        next_ele(1) = la_ele(1);
        next_ele(2) = la_ele(3);
        next_azi(1) = 0;
        next_azi(2) = 2*pi - 2*pi/n_rays_azi;
    end

    la_ele = linspace(next_ele(1), next_ele(2), n_rays_ele);
    la_azi = linspace(next_azi(1), next_azi(2), n_rays_azi);

    ii = ii+1;

    % pause(0.01)
end

hit = best_d <= rad_rcv;
best_r = r(:,:,best_ray);

end

