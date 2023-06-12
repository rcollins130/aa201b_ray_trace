function [best_r] = find_eigenvector(dt, nt, dz, x_src,x_rcv, rad_rcv, z,v,c)
%FIND_EIGENVECTOR Summary of this function goes here
%   Detailed explanation goes here

% initial sweep pattern
n_rays_ele_init = 7;
n_rays_azi_init = 8;

n_rays_ele = 5;
n_rays_azi = 5;

la_ele = linspace(-pi/2, pi/2, n_rays_ele_init);
la_azi = linspace(0, 2*pi - 2*pi/n_rays_azi, n_rays_azi_init);

best_d = 9e9;

while best_d > rad_rcv
    
    % construct rays (replicates launch_angles)
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
    
    % trace rays
    r = zeros(nt, 3, size(la, 1));
    max_ii_t = repmat(nt, size(la,1),1);
    
    for ii_ray = 1:size(la,1)
        [r(:,:,ii_ray), max_ii_t(ii_ray)]= raytrace_rac( ...
            dt, nt, dz, x_src, la(ii_ray,1), la(ii_ray,2), z, v, c ...
            );
    end
    
    % find ray that got closest to the point
    d_mins = min(vecnorm(x_rcv - r, 2,2),[],1);
    [best_d, best_ray] = min(d_mins)
    
    % get neighbors of best ray
    next_ele_idx = mod(ele_idx(best_ray) + [-1,1] -1, length(la_ele))+1;
    next_azi_idx = mod(azi_idx(best_ray) + [-1,1] -1, length(la_azi))+1;
    next_ele = la_ele(next_ele_idx);
    next_azi = la_azi(next_azi_idx);
    if next_azi(1) > next_azi(2)
        next_azi(1) = next_azi(1) - 2*pi;
    end

    la_ele = linspace(next_ele(1), next_ele(2), n_rays_ele);
    la_azi = linspace(next_azi(1), next_azi(2), n_rays_azi);
end

best_r = r(:,:,best_ray);

end

