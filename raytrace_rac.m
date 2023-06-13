function [this_r, ii_t] = raytrace_rac(dt, nt, dz, r_0, la_ele, la_azi, z, v, c)
%RAYTRACE_RAC Summary of this function goes here
%   INPUT
%   dt: time step, s
%   r_0: [x,y,z] Initial value
%   la_ele: Launch Angle Elevation, from XY
%   la_azi: Launch Angle Azimuth, from XZ
%   z: [z;] Altitude Profile
%   v: [vx, vy] Velocity Profile
%   c: [c;] Sound Speed Profile
%   OUTPUT
%   r: [x,y,z;] Ray Position
%   ii_t: time index of last valid point on ray

% compute gradients
[~, v_dz] = gradient(v, dz);
c_dz = gradient(c, dz);

% Initialize ray position vector
this_r = nan(nt, 3);
this_r(1,:) = r_0;

% Initialize Ray Normal
n = [
    cos(la_ele) * cos(la_azi) ...
    -cos(la_ele) * sin(la_azi) ...
    sin(la_ele)
];

% Function for determining profile index from z value
get_zidx = @(zt)(floor((abs(zt - z(1))) / dz)+1);

% Get initial z index
ii_z1 = get_zidx(r_0(3)); 

% compute initial slowness vector
s = n ./ (c(ii_z1) + n * [v(ii_z1,:), 0]');

% compute slowness perp and omega vec
s_perp = s(1:2);
s_z = s(3);
omega = 1 - v(:,1:2) * s(1:2)';

% setup RK4 function
f = @(sgn, ii_z, s_z) [
    c(ii_z).^2 ./ omega(ii_z) * s_perp(1) + v(ii_z,1); 
    c(ii_z).^2 ./ omega(ii_z) * s_perp(2) + v(ii_z,2); 
    c(ii_z).^2 ./ omega(ii_z) * s_z;
    sgn*(-omega(ii_z)./c(ii_z) * c_dz(ii_z) - s_perp * v_dz(ii_z,1:2)');
];

% iterate thru time
for ii_t = 1:size(this_r,1)-1
    % break if beyond z bounds or time reached
    if ii_z1 >= size(z, 1) || ii_t == size(this_r,1)
        break
    end
    
    % compute rk4 coefficents
    z1 = this_r(ii_t, 3);
    a1 = f(sign(z1), ii_z1, s_z);

    z2 = z1 + dt/2 * a1(3);
    ii_z2 = get_zidx(z2);
    if ii_z2 > size(z,1); ii_z2 = size(z,1); end
    a2 = f(sign(z2), ii_z2, s_z + dt/2 * a1(4));

    z3 = z1 + dt/2 * a2(3);
    ii_z3 = get_zidx(z3);
    if ii_z3 > size(z,1); ii_z3 = size(z,1); end

    a3 = f(sign(z3), ii_z3, s_z + dt/2 * a2(4));

    z4 = z1 + dt * a3(3);
    ii_z4 = get_zidx(z4);
    if ii_z4 > size(z,1); ii_z4 = size(z,1); end
    a4 = f(sign(z4), ii_z4, s_z + dt * a3(4));
    
    % compute dx
    dx = dt/6 * (a1 + 2*a2 + 2*a3 + a4);
    
    % advance r and s_z
    this_r(ii_t+1, :) = this_r(ii_t, :) + dx(1:3)';
    s_z = s_z + dx(4);
    
    % get new z index
    ii_z1 = get_zidx(this_r(ii_t+1,3));
end

end
