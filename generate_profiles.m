function [T,c,v] = generate_profiles(z, wind_azi)
%GENERATE_PROFILES Summary of this function goes here
%   Detailed explanation goes here

T_0 = 288.15; % Reference sea level temperature, K
L_b = 0.0065; % Temperature lapse rate, K/m
gamma= 1.4; % adiabatic index
R = 287.058; % molar gas constant, dry air, m2/s2
kar = 0.40; % Karman Constant
z_0 = 0.1; % Surface roughness, m 
v_s = 0.6; % Friction velocity, m/s

% temperature profile
T = T_0 - L_b * z;

% sound speed profile
c = sqrt(gamma * R * T);

if nargin > 1
% velocity profile
% TODO: ADD POWER LAW ABOVE ~1000M?
v_mag= v_s/kar .* log(z/z_0);
v = v_mag * [cos(wind_azi), sin(wind_azi)]; % x, y
v(z<z_0, :) = 0;
end

end

