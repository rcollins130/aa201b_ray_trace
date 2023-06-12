

% demo ray tracing and finding eigenvectors from moving source to 
% observer

% set up environemnt
dt = 0.1; %small time, s
t_max = 20; 
nt = t_max / dt+1;

dz = 0.1;
z_range=[0,2000];

wind_azi = deg2rad(0);
z = linspace(z_range(1), z_range(2), (diff(z_range))/dz + 1)';
[T,c,v] = generate_profiles(z, wind_azi);
v = v;

% observer point
x_rcv = [50,1000,10];

% aircraft path, from paper
act = ArntzenAircraftTrajectory(1);

% plot scene
figure(2); clf; hold on;
plot3(x_rcv(1), x_rcv(2), x_rcv(3), 'bo', 'MarkerFaceColor','b');
plot3(act.position(:,1), act.position(:,2), act.position(:,3))
xlabel('x')
ylabel('y')
zlabel('z')

xlim([-4000, 4000]);
ylim([-4000, 4000]);
zlim([0, 1000]);

ap = plot3(act.position(1,1), act.position(1,2), act.position(1,3),'ro','MarkerFaceColor','r');

% loop thru each aircraft time 
for ii_bigt = 1:5:size(act.time, 1)
    % get eigenray from point to 
   [egv, hit] = find_eigenvector(dt, nt, dz, act.position(ii_bigt,:),x_rcv, 15, z, v, c);
    
    if hit
        plot3(egv(:,1), egv(:,2), abs(egv(:,3)),'g');
    else
        plot3(egv(:,1), egv(:,2), abs(egv(:,3)),'g--');
    end
    ap.XData = act.position(ii_bigt, 1);
    ap.YData = act.position(ii_bigt, 2);
    ap.ZData = act.position(ii_bigt, 3);

    pause(0.01);
end