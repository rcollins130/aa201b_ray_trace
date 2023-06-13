%%
base_path ="/Users/robertcollins/Stanford/Stanford_Google_Drive/AA201B/Final_Project/";
dir = "figures";
outpath = fullfile(base_path,dir);
mkdir(outpath);

%% Environment Figures

dz = 1; %m
z_range = [0, 3500]; %m
wind_azi = 0;

z = linspace(z_range(1), z_range(2), (diff(z_range))/dz + 1)';
[T,c,v] = generate_profiles(z, wind_azi);

% 
figure(1); clf;
subplot(1,3,1)
plot(T,z)
xlabel('Temperature, K')
ylabel('Altitude, m')

subplot(1,3,2)
plot(c,z)
xlabel('Sound Speed, m/s')
ylabel('Altitude, m')

subplot(1,3,3)
plot(vecnorm(v,2,2),z)
xlabel('Wind Velocity, m/s')
ylabel('Altitude, m')

sgtitle('Standard Profile')
saveas(gcf, fullfile(outpath, 'profile_standard.png'))

% 
c = -sawtooth(2*pi*z/max(z), 1/2)*10 + 350;

figure(1); clf;
subplot(1,3,1)
plot(T,z)
xlabel('Temperature, K')
ylabel('Altitude, m')

subplot(1,3,2)
plot(c,z)
xlabel('Sound Speed, m/s')
ylabel('Altitude, m')

subplot(1,3,3)
plot(vecnorm(v,2,2),z)
xlabel('Wind Velocity, m/s')
ylabel('Altitude, m')

sgtitle('Duct Profile')
saveas(gcf, fullfile(outpath, 'profile_duct.png'))

% static 
v = v*0 ;
c = 343*ones(size(z));

figure(1); clf;
subplot(1,3,1)
plot(T,z)
xlabel('Temperature, K')
ylabel('Altitude, m')

subplot(1,3,2)
plot(c,z)
xlabel('Sound Speed, m/s')
ylabel('Altitude, m')

subplot(1,3,3)
plot(vecnorm(v,2,2),z)
xlabel('Wind Velocity, m/s')
ylabel('Altitude, m')

sgtitle('Static Profile')
saveas(gcf, fullfile(outpath, 'profile_static.png'))