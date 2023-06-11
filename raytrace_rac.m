function [this_r, ii_t] = raytrace_rac(dt, nt, dz, r_0, la_ele, la_azi, z, v, c)
%RAYTRACE_RAC Summary of this function goes here
%   INPUT
%   dt: time step, s
%   r_0: [x,y,z] Initial 
%   la_ele: Launch Angle Elevation, from XY
%   la_azi: Launch Angle Azimuth, from XZ
%   z: [z;] Altitude Profile
%   v: [vx, vy] Velocity Profile
%   c: [c;] Sound Speed Profile
%   OUTPUT
%   r: [x,y,z;] Ray Position

% compute gradients
v_dz = [diff(v,1,1)/dz; zeros(1,2)];
c_dz = [diff(c,1,1)/dz; 0];

% Initialize ray position vector
this_r = zeros(nt, 3);
this_r(1,:) = r_0;

% Initialize Ray Normal
n = [
    cos(la_ele) * cos(la_azi) ...
    -cos(la_ele) * sin(la_azi) ...
    sin(la_ele)
];

% compute starting ii_z
ii_z0 = round((abs(r_0(3) - z(1))) / dz); 

% [~, ii_z0] = min(abs(z-this_r(1, 3)));

% compute initial slowness vector
s = n ./ (c(ii_z0) + n * [v(ii_z0,:), 0]');

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

ii_z1 = ii_z0;
% iterate thru time
for ii_t = 1:size(this_r,1)-1
    if ii_z1 >= size(z, 1) || ii_t == size(this_r,1) %|| ...
        % this_r(ii_t,1) > x_range(2) || ...
        % this_r(ii_t,1) < x_range(1) || ...
        % this_r(ii_t,2) > y_range(2) || ...
        % this_r(ii_t,2) < y_range(1)
        break
    end
    
    % if use_rk4
        % advance values (RK4)
        z1 = this_r(ii_t, 3);
        a1 = f(sign(z1), ii_z1, s_z);

        z2 = z1 + dt/2 * a1(3);
        ii_z2 = abs(ii_z1 + ceil((abs(z2)-z1)/dz));
        if ii_z2 > size(z,1); break; end
        a2 = f(sign(z2), ii_z2, s_z + dt/2 * a1(4));

        z3 = z1 + dt/2 * a2(3);
        ii_z3 = abs(ii_z1 + ceil((abs(z3)-z1)/dz));
        if ii_z3 > size(z,1); break; end

        a3 = f(sign(z3), ii_z3, s_z + dt/2 * a2(4));

        z4 = z1 + dt * a3(3);
        ii_z4 = abs(ii_z1 + ceil((abs(z4)-z1)/dz));
        if ii_z4 > size(z,1); break; end
        a4 = f(sign(z4), ii_z4, s_z + dt * a3(4));

        dx = dt/6 * (a1 + 2*a2 + 2*a3 + a4);

        this_r(ii_t+1, :) = this_r(ii_t, :) + dx(1:3)';
        s_z = s_z + dx(4);
    % else
    %     % advance values (euler method)
    %     warning('reflections not implemented');
    %     this_r(ii_t+1, :) = this_r(ii_t,:) + ...
    %         (c(ii_z1).^2 / omega(ii_z1) .* [s_perp, s_z] + v(ii_z1,:)) .* dt;
    %     s_z = s_z + sz_dt(ii_z1) * dt;
    % end
    
    ii_z1 = ii_z0 + ceil((abs(this_r(ii_t+1, 3))-r_0(3))/dz);
end

end
