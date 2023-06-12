function [la, ele_idx, azi_idx] = launch_angles(min_ele,max_ele, n_rays_ele, min_azi,max_azi,n_rays_azi)
%LAUNCH_ANGLES Summary of this function goes here
%   Detailed explanation goes here

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

end

