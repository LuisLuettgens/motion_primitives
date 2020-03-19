clear
% close
figure;
hold on

% [A,map] = imread('occ_map_wo_buffer.ppm');
% A_compl = imcomplement(A);
% lat=[54.0906409 54.235];
% lon=[12.0416667 12.1547041];
lat=[-8895.613288 7172.437734];
lon=[-3865.449476 3518.68377];
% image([lon(1)+0.5*(lon(2)-lon(1))/size(A,1) lon(2)-0.5*(lon(2)-lon(1))/size(A,1)],...
%       [lat(1)+0.5*(lat(2)-lat(1))/size(A,2) lat(2)-0.5*(lat(2)-lat(1))/size(A,2)],permute(A_compl,[2 1 3]));
  
%   Latitude: (lat(end)-lat(1))*110.574
% Longitude: (lon(end)-lon(1))*111.320*cos( (lat(1)+lat(end))/2. )
  
axis equal;
% set(gca,'DataAspectRatio',[110.574 111.320*cosd((lat(1)+lat(end))/2.) 1]);

lat0 = 54.17057475;
lon0 = 12.10074142;
harbor_segs=textread('RostockHarbor/harbor.txt');
harbor_segs_ned=harbor_segs;
for i=1:length(harbor_segs)
    [harbor_segs_ned(i,1),harbor_segs_ned(i,2),~]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0.,lon0,lat0,0.);
end
plot(harbor_segs_ned(:,2),harbor_segs_ned(:,1),'k', 'LineWidth',1.5);

harbor_segs=textread('RostockHarbor/obstacle_01.txt');
harbor_segs_ned=harbor_segs;
for i=1:length(harbor_segs)
    [harbor_segs_ned(i,1),harbor_segs_ned(i,2),~]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0.,lon0,lat0,0.);
end
plot(harbor_segs_ned(:,2),harbor_segs_ned(:,1),'k', 'LineWidth',1.5);

harbor_segs=textread('RostockHarbor/obstacle_02.txt');
harbor_segs_ned=harbor_segs;
for i=1:length(harbor_segs)
    [harbor_segs_ned(i,1),harbor_segs_ned(i,2),~]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0.,lon0,lat0,0.);
end
plot(harbor_segs_ned(:,2),harbor_segs_ned(:,1),'k', 'LineWidth',1.5);

harbor_segs=textread('RostockHarbor/obstacle_03.txt');
harbor_segs_ned=harbor_segs;
for i=1:length(harbor_segs)
    [harbor_segs_ned(i,1),harbor_segs_ned(i,2),~]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0.,lon0,lat0,0.);
end
plot(harbor_segs_ned(:,2),harbor_segs_ned(:,1),'k', 'LineWidth',1.5);

harbor_segs=textread('RostockHarbor/obstacle_04.txt');
harbor_segs_ned=harbor_segs;
for i=1:length(harbor_segs)
    [harbor_segs_ned(i,1),harbor_segs_ned(i,2),~]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0.,lon0,lat0,0.);
end
plot(harbor_segs_ned(:,2),harbor_segs_ned(:,1),'k', 'LineWidth',1.5);

harbor_segs=textread('RostockHarbor/obstacle_05.txt');
harbor_segs_ned=harbor_segs;
for i=1:length(harbor_segs)
    [harbor_segs_ned(i,1),harbor_segs_ned(i,2),~]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0.,lon0,lat0,0.);
end
plot(harbor_segs_ned(:,2),harbor_segs_ned(:,1),'k', 'LineWidth',1.5);

harbor_segs=textread('RostockHarbor/obstacle_06.txt');
harbor_segs_ned=harbor_segs;
for i=1:length(harbor_segs)
    [harbor_segs_ned(i,1),harbor_segs_ned(i,2),~]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0.,lon0,lat0,0.);
end
plot(harbor_segs_ned(:,2),harbor_segs_ned(:,1),'k', 'LineWidth',1.5);

harbor_segs=textread('RostockHarbor/obstacle_07.txt');
harbor_segs_ned=harbor_segs;
for i=1:length(harbor_segs)
    [harbor_segs_ned(i,1),harbor_segs_ned(i,2),~]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0.,lon0,lat0,0.);
end
plot(harbor_segs_ned(:,2),harbor_segs_ned(:,1),'k', 'LineWidth',1.5);

harbor_segs=textread('RostockHarbor/obstacle_08.txt');
harbor_segs_ned=harbor_segs;
for i=1:length(harbor_segs)
    [harbor_segs_ned(i,1),harbor_segs_ned(i,2),~]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0.,lon0,lat0,0.);
end
plot(harbor_segs_ned(:,2),harbor_segs_ned(:,1),'k', 'LineWidth',1.5);

harbor_segs=textread('RostockHarbor/obstacle_09.txt');
harbor_segs_ned=harbor_segs;
for i=1:length(harbor_segs)
    [harbor_segs_ned(i,1),harbor_segs_ned(i,2),~]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0.,lon0,lat0,0.);
end
plot(harbor_segs_ned(:,2),harbor_segs_ned(:,1),'k', 'LineWidth',1.5);

harbor_segs=textread('RostockHarbor/obstacle_10.txt');
harbor_segs_ned=harbor_segs;
for i=1:length(harbor_segs)
    [harbor_segs_ned(i,1),harbor_segs_ned(i,2),~]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0.,lon0,lat0,0.);
end
plot(harbor_segs_ned(:,2),harbor_segs_ned(:,1),'k', 'LineWidth',1.5);

run visited_positions.m
plot(vis_pos(:,2),vis_pos(:,1),'b.')

% run points.m
% plot(pts(:,2),pts(:,1),'b.')

% for i=1:length(vis_pos)
%     theta = vis_pos(i,1);
%     r = 10; % magnitude (length) of arrow to plot
%     x = vis_pos(i,2); y = vis_pos(i,3);
%     u = r * cos(theta); % convert polar (theta,r) to cartesian
%     v = r * sin(theta);
%     h = quiver(y,x,v,u,'Color',[0,0,1],'MarkerSize',1,'LineWidth',1,'MaxHeadSize',1);
% end

% run route_points.m
% plot(route(:,1),route(:,2),'r.')

run solution.m
plot(sol(:,2),sol(:,1),'r-','Markersize',8, 'LineWidth',2.5)
plot(sol(1,2),sol(1,1),'ro','Markersize',10);
s = text(sol(1,2),sol(1,1)+100,'Start');
%s.FontWeight = 'bold'
s.FontSize = 14;
plot(sol(end,2),sol(end,1),'ro','Markersize',10)
t = text(sol(end,2)+130,sol(end,1)+125,'Ziel');
%t.FontWeight = 'bold'
t.FontSize = 14;
%a = load("intersection_points.m")
%plot(a(:,2),a(:,1),'b.')
%plot(intersection_points(:,2),intersection_points(:,1),'b.')
axis([-3000 2000 -2500 3750])

%xlabel('$y$ [m]','Interpreter','latex')
%ylabel('$x$ [m]','Interpreter','latex')

set(gcf,'Position',[356          31        1280         873]);

print('arc_traj','-depsc','-r600');
% for i=1:length(route)
%     theta = route(i,1);
%     r = 10; % magnitude (length) of arrow to plot
%     x = route(i,2); y = route(i,3);
%     u = r * cos(theta); % convert polar (theta,r) to cartesian
%     v = r * sin(theta);
%     h = quiver(y,x,v,u,'Color',[0,0,0],'MarkerSize',1,'LineWidth',1,'MaxHeadSize',1);
%     
%     u_ = r * cos(theta+alpha); % convert polar (theta,r) to cartesian
%     v_ = r * sin(theta+alpha);
%     h = quiver(y,x,v_,u_,'Color',[0,1,0],'MarkerSize',4,'LineWidth',2,'MaxHeadSize',3);
% end