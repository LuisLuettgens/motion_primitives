clear
close

figure;
hold on
axis equal

long0=12.10074142;
lat0=54.17057475;
h0=0;

harbor_segs=textread('harbor_polygon_71vertices.txt')

for i=1:length(harbor_segs)
    [x,y,z]=convertLLAtoNED(harbor_segs(i,1),harbor_segs(i,2),0,long0,lat0,h0);
    harbor_segs(i,1)=x;
    harbor_segs(i,2)=y;
end

%harbor_segs = csvread('poly_simpler_xy.csv');

max_x=-inf;
max_y=-inf;
min_x= inf;
min_y= inf;
for i=1:length(harbor_segs)
    if harbor_segs(i,1)>max_x
        max_x=harbor_segs(i,1);
    end
    if harbor_segs(i,2)>max_y
        max_y=harbor_segs(i,2);
    end
    
    if harbor_segs(i,1)<min_x
        min_x=harbor_segs(i,1);
    end
    if harbor_segs(i,2)<min_y
        min_y=harbor_segs(i,2);
    end
end


%grid.width=310;
grid.width=0.003;

min_x=min_x-grid.width/8.;
min_y=min_y-grid.width/8.;

plot(harbor_segs(:,1),harbor_segs(:,2));

grid.x=[];
grid.y=[];

x=min_x;
grid.min_x=x;
grid.x=[grid.x x];
grid.x_size=0;
while x<=max_x
    x=x+grid.width;
    grid.x=[grid.x x];
    grid.x_size=grid.x_size+1;
end
grid.max_x=x;
y=min_y;
grid.min_y=y;
grid.y=[grid.y y];
grid.y_size=0;
while y<=max_y
    y=y+grid.width;
    grid.y=[grid.y y];
    grid.y_size=grid.y_size+1;
end
grid.max_y=y;

x=min_x;
plot([x x],[min_y grid.max_y],'k');
while x<=max_x
    x=x+grid.width;
    plot([x x],[min_y grid.max_y],'k');
end

y=min_y;
plot([min_x grid.max_x],[y y],'k');
while y<=max_y
    y=y+grid.width;
    plot([min_x grid.max_x],[y y],'k');
end

m_pt1=[0 0];
head1=7.*pi/4.;
leng1=300.;
wid1=75.;

% seg=get_ship_segs(m_pt1,head1,leng1,wid1);

for i=1:grid.x_size
    fprintf('calculating: %2i of %2i\n',i,grid.x_size);
    for j=1:grid.y_size
        cell_segs=[grid.x(i) grid.y(j); grid.x(i+1) grid.y(j); grid.x(i+1) grid.y(j+1); grid.x(i) grid.y(j+1)];
        [A, clippedPolygon]=calc_clipped_area(harbor_segs,cell_segs);
        if clippedPolygon
            grid.harbor_area{i,j}.x=clippedPolygon(:,1)';
            grid.harbor_area{i,j}.y=clippedPolygon(:,2)';
            if mod(i+j,2)
                p2=patch(clippedPolygon(:,1),clippedPolygon(:,2),'g');
                set(p2,'EdgeAlpha',0);
            else
                p2=patch(clippedPolygon(:,1),clippedPolygon(:,2),'r');
                set(p2,'EdgeAlpha',0);
            end
        end
        
        grid.segs_in_dist{i,j}=[];
        
        nearest_seg=0;
        nearest_dist=inf;
        
        for k=1:length(harbor_segs)-1

            for l=1:length(cell_segs)-1

                [dist_,s,t,c1,c2]=ClosestPtSegmentSegment(cell_segs(l,:), cell_segs(l+1,:), harbor_segs(k,:), harbor_segs(k+1,:));

                if dist_<nearest_dist
                    nearest_dist=dist_;
                    nearest_seg=k;
                end
                
            end
                
        end
        
        max_dist_of_nearest_seg=-inf;
        
        for k=nearest_seg:nearest_seg+1

            for l=1:length(cell_segs)
                
                dist_=(harbor_segs(k,1)-cell_segs(l,1))^2 + (harbor_segs(k,2)-cell_segs(l,2))^2;
                
                if dist_>max_dist_of_nearest_seg
                    
                    max_dist_of_nearest_seg=dist_;
                    
                end
                
            end
            
        end
        
        for k=1:length(harbor_segs)-1

            for l=1:length(cell_segs)-1

                [dist_,s,t,c1,c2]=ClosestPtSegmentSegment(cell_segs(l,:), cell_segs(l+1,:), harbor_segs(k,:), harbor_segs(k+1,:));

                if dist_<max_dist_of_nearest_seg
                    grid.segs_in_dist{i,j}=[grid.segs_in_dist{i,j} k];
                    break;
                end

            end

        end
        
    end
end

xml_write('grid_71vertices.xml', grid);

