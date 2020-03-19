clear
close

level=-8;
center_of_harbor=[12.1034 54.16];

cloud1 = readtable('bathym_data/EPSG_4326-650-6000-50.txt');
cloud2 = readtable('bathym_data/EPSG_4326-700-5950-50.txt');
cloud3 = readtable('bathym_data/EPSG_4326-700-6000-50.txt');

harb_f={'input_polygons/harbor_polygon_w_o_obst_lla.txt'};
obst_f={'input_polygons/obst1.txt';...
        'input_polygons/obst2.txt'};

harbor_segs=textread(harb_f{1});


x=[cloud1.x; cloud2.x; cloud3.x];
y=[cloud1.y; cloud2.y; cloud3.y];
z=[cloud1.z; cloud2.z; cloud3.z];

sc_int = scatteredInterpolant(x,y,z,'linear','none');
x_arr=linspace(min(x),max(x),400);
y_arr=linspace(min(y),max(y),800);

[qx,qy] = meshgrid(x_arr,y_arr);
qz = sc_int(qx,qy);
for i=1:size(qz,1)
    for j=1:size(qz,2)
        if isnan(qz(i,j))
            qz(i,j)=0;
        end
    end
end

% for i=1:size(qx,1)
%     for j=1:size(qx,2)
%         if ~pnpoly(harbor_segs(:,1),harbor_segs(:,2),qx(i,j),qy(i,j))
%             if isnan(qz(i,j))
%                 qz(i,j)=-20.;
%             elseif qz(i,j)<=level
%                 qz(i,j)=-20.;
%             else
%                 qz(i,j)=0.;
%             end
%         elseif isnan(qz(i,j))
%             qz(i,j)=-20.;
%         end
%     end
% end

% gr_int = griddedInterpolant(qx',qy',qz');
% 
% mean_err=0;
% max_err=-inf;
% for i=1:length(x)
%     act_err=abs(gr_int(x(i),y(i))-z(i));
%     if ~isnan(act_err)
%         mean_err=mean_err+act_err;
%         if act_err>max_err
%             max_err=act_err;
%         end
%     end
% end
% mean_err=mean_err/length(x);
% 
% mean_err
% max_err

figure;
hold on
% axis equal

% plot3(harbor_segs(:,1),harbor_segs(:,2),zeros(length(harbor_segs),1),'k');
plot(harbor_segs(:,1),harbor_segs(:,2),'k');
plot(center_of_harbor(1),center_of_harbor(2),'go','MarkerSize',15)
plot(center_of_harbor(1),center_of_harbor(2),'gx','MarkerSize',15)

% scatter3(x,y,z,[],zeros(length(x),1),'.')

set(gca, 'DataAspectRatio', [repmat(max(diff(get(gca, 'XLim')), diff(get(gca, 'YLim'))), [1 2]) diff(get(gca, 'ZLim'))])

% view([0,0])

% figure;
% mesh(qx,qy,qz);

[C,h]=contour(qx,qy,qz,[level level]);%,'ShowText','on');

set(gca,'XLim',[x_arr(1)-abs(x_arr(end)-x_arr(1))*0.05 x_arr(end)+abs(x_arr(end)-x_arr(1))*0.05])
set(gca,'YLim',[y_arr(1)-abs(y_arr(end)-y_arr(1))*0.05 y_arr(end)+abs(y_arr(end)-y_arr(1))*0.05])

% save temporary contour file
fileID = fopen('contour_lines_temp.txt','w');

i=1;
while i<=length(C)
    n=C(2,i);
    j=i+1;
    while j<=i+n
        fprintf(fileID,'%12.8f ',C(1,j));
        fprintf(fileID,'%12.8f\n',C(2,j));
        j=j+1;
    end
    i=j;
    fprintf(fileID,'\n');
end
fclose(fileID);

status = system(['./polyintersec ' harb_f{1} ' '...
                 'contour_lines_temp.txt '...
                 'intersec_poly_temp.txt']);
                        
status = system('rm contour_lines_temp.txt');
                        
fileID = fopen(['intersec_poly_temp.txt'],'r');
tline = fgets(fileID);
i=1;
count=1;
while ischar(tline)
    if strcmp(tline,char(10))
        count=count+1;
        i=1;
    else
        deep_water{count}(i,:)=strread(tline, '%f', 2)';
        i=i+1;
    end
    tline = fgets(fileID);
end
fclose(fileID);

status = system('rm intersec_poly_temp.txt');

% the variable "center_of_harbor" defines the main harbor-area:
% go over all "deep_water"-polygons to find the one with "center_of_harbor"
% in it
found=false;
harbor_poly_nr=0;
i=1;
while i<=length(deep_water) && ~found
    if pnpoly(deep_water{i}(:,1),deep_water{i}(:,2),center_of_harbor(1),center_of_harbor(2))
        found=true;
        harbor_poly_nr=i;
    end
    i=i+1;
end
plot([deep_water{harbor_poly_nr}(:,1);deep_water{harbor_poly_nr}(1,1)],[deep_water{harbor_poly_nr}(:,2);deep_water{harbor_poly_nr}(1,2)],'r--');

% save harbor file
fileID = fopen(['output_polygons/harbor_deeper_than_' num2str(level) '.txt'],'w');
for j=1:length(deep_water{harbor_poly_nr})
   fprintf(fileID,'%12.8f ',deep_water{harbor_poly_nr}(j,1));
   fprintf(fileID,'%12.8f\n',deep_water{harbor_poly_nr}(j,2));
end
fclose(fileID);

% save all obstacles which are inside the harbor deep water
obstacle_count=0;
for i=[1:harbor_poly_nr-1 harbor_poly_nr+1:length(deep_water)]
    % its enough to check one point of each obstacle, because the polygons
    % are not intersecting
    if pnpoly(deep_water{harbor_poly_nr}(:,1),deep_water{harbor_poly_nr}(:,2),deep_water{i}(1,1),deep_water{i}(1,2))
        obstacle_count=obstacle_count+1;
        plot([deep_water{i}(:,1);deep_water{i}(1,1)],[deep_water{i}(:,2);deep_water{i}(1,2)],'b--');
        %save obstacle file
        fileID = fopen(['output_polygons/obstacle_shallower_than_' num2str(level) '_nr_' num2str(obstacle_count) '.txt'],'w');
        for j=1:length(deep_water{i})
            fprintf(fileID,'%12.8f ',deep_water{i}(j,1));
            fprintf(fileID,'%12.8f\n',deep_water{i}(j,2));
        end
        fclose(fileID);
    end
end