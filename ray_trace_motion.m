

% demo ray tracing and finding eigenvectors from moving source to 
% observer

% set up environemnt
dt = 0.01; %small time, s
t_max = 20; 
nt = t_max / dt+1;

dz = 1;
z_range=[0,1000];

rot=1;
wad = 90;
wind_azi = deg2rad(wad);
z = linspace(z_range(1), z_range(2), (diff(z_range))/dz + 1)';
[T,c,v] = generate_profiles(z, wind_azi);
v = 0*v;

% observer point
x_rcv = [50,2000,2];

% aircraft path, from paper
traj = 2;
act = ArntzenAircraftTrajectory(traj);

% plot scene
figure(2); clf; hold on;
box on; 
grid on;
view([-1,1,1])
plot3(x_rcv(1), x_rcv(2), x_rcv(3), 'bo', 'MarkerFaceColor','b');
plot3(act.position(:,1), act.position(:,2), act.position(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
title(sprintf("Trajectory %d, Wind Azimuth %0.1f", traj, wad))
% title(sprintf("Trajectory %d, no wind", traj))


xlim([-4000, 4000]);
ylim([-4000, 4000]);
zlim([0, 1000]);

ap = plot3(act.position(1,1), act.position(1,2), act.position(1,3),'ro','MarkerFaceColor','r');
er = plot3(zeros(nt, 1), zeros(nt, 1), zeros(nt, 1), 'b','LineWidth',2);

ii=1;
% loop thru each aircraft time 
for ii_bigt = 1:5:size(act.time, 1)
    % get eigenray from point to 
   [egv, hit] = find_eigenvector(dt, nt, dz, act.position(ii_bigt,:),x_rcv, 5, z, v, c);
    
    er.XData = egv(:,1);
    er.YData = egv(:,2);
    er.ZData = abs(egv(:,3));
    if hit 
        er.LineStyle = '-';
    else
        er.LineStyle = '--';
    end
    
    msk = egv(:,3)>0;
    if hit
        plot3(egv(msk,1), egv(msk,2), abs(egv(msk,3)),'g');
    else
        plot3(egv(msk,1), egv(msk,2), abs(egv(msk,3)),'g--');
    end
    ap.XData = act.position(ii_bigt, 1);
    ap.YData = act.position(ii_bigt, 2);
    ap.ZData = act.position(ii_bigt, 3);

    pause(0.01);
    % [az,el] = view;
    % view(az+5,el);
    f(ii) = getframe(gcf);
    ii=ii+1;
end

if rot
    for ii_t = 1:360
        [az,el] = view;
        view(az+1,el);
        pause(0.1)
        f(ii) = getframe(gcf);
        ii=ii+1;
    end
end

v = VideoWriter(sprintf('vid_plane_%d_%d_%d',traj,wad,rot),'MPEG-4');
open(v)
writeVideo(v,f)
close(v)
clear f