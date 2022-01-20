global obstacles
global world
global repulsiveWeight
global d_star;
global zeta;
global start;
global steps;
global epsilon;
global goal;

d_star = 10;
start = [-5,-15,-10,0,-5,-5;-4,5,12,3,-15,-12];
steps = 1000;
goal = [15;0];
repulsiveWeight = 2;
zeta = 0.1;
epsilon = 0.1;
world = struct('center',[0;0], 'radius', 20,'distanceInfluence',2);
obstacles = [struct('center',[-8;5],'radius',5,'distanceInfluence',4); struct('center',[0;-9],'radius',4,'distanceInfluence',8); struct('center',[5;0],'radius',4,'distanceInfluence',8)];

close all;
plot_surface()
plot_contour()

%% Plotting
function plot_contour()
global obstacles;
global goal;
global start;

[xx, yy] = generate_meshgrid();
[UAttr,URep] = calculate_total_potential(xx,yy);

maximum = max(UAttr, [], 'all');
UTotal = UAttr + URep;
UTotal(UTotal > maximum) = maximum;

figure
contour(xx, yy, UTotal, 50)
hold on

plot(goal(1),goal(2),'r*');

for i = 1:numel(obstacles)
    obstacle = obstacles(i);
    plot(obstacle.center(1,1), obstacle.center(2,1));
    hold on
end

if size(start,2) > 0
        colors = ['r', 'b', 'g', 'm', 'c','k'];
    for i = 1:size(start,2)
        start_ = start(:,i);
        path = potential_planner(start_);
        color = "-"+colors(i);
        plot(start_(1), start_(2),"*"+color);
        plot(path(1,:), path(2,:),"."+color);
    end
else
    path = potential_planner(start);
    plot(path(1,:), path(2,:), '-r')
end
end

function plot_surface()
figure
x = linspace(-25, 25);
y =linspace(-25, 25);
[xx, yy] = meshgrid(x, y);

[UAttr,URep] = calculate_total_potential(xx, yy);


maximum = max(UAttr, [], 'all');
UTotal =  URep + UAttr ;
UTotal(UTotal > maximum) = maximum;
surf(xx, yy, UTotal,'FaceColor','interp','EdgeColor','interp')
set(gca,'XTick',[],'YTick',[],'ZTick',[])
title('Total Potential Surface')

figure
UTotal = UAttr ;
UTotal(UTotal > maximum) = maximum;
surf(xx, yy, UTotal,'FaceColor','interp','EdgeColor','interp')
set(gca,'XTick',[],'YTick',[],'ZTick',[])
title('Attractive Potential Surface')

figure
UTotal = URep ;
UTotal(UTotal > maximum) = maximum;
surf(xx, yy, UTotal,'FaceColor','interp','EdgeColor','interp')
set(gca,'XTick',[],'YTick',[],'ZTick',[])
title('Repulsive Potential Surface')
end

function [UAttr,URep] = calculate_total_potential(xx,yy)
UAttr = [];
URep = [];

for i = 1:numel(xx)
    qx = xx(i);
    qy = yy(i);
    q = [qx; qy];
    attractive = potential_attractive(q);
    repulsive = potential_repulsiveSphere(q);
    UAttr(i) = attractive;
    URep(i) = repulsive;
end

UAttr = reshape(UAttr, size(xx));
URep = reshape(URep, size(xx));

end

function [xx,yy] = generate_meshgrid()
global world;

x_min = world.center(1) - world.radius;
x_max = world.center(1) + world.radius;
y_min = world.center(2) - world.radius;
y_max = world.center(2) + world.radius;
x = linspace(x_min - 0.1,x_max + 0.1);
y = linspace(y_min - 0.1, y_max + 0.1);

[X,Y] = meshgrid(x, y);
xx = X;
yy = Y;
end


%% Repulsive Potential
function [URep] = potential_repulsiveSphere(q)
global world;
global repulsiveWeight;
global obstacles;

URep = 0;
dist = world.radius - distance(q, world.center);
if dist < 0.0005
    URep = 1/0;
elseif abs(dist) <= world.distanceInfluence
    URep = repulsiveWeight * (1/2) * ((1/dist) - (1/world.distanceInfluence))^2;
end

for i = 1:numel(obstacles)
    obstacle = obstacles(i);
    obstacle_distance = distance_toSphere(q, obstacle);
    if obstacle_distance < 0.0005
        URep = URep + (1/0);
    elseif obstacle_distance <= obstacle.distanceInfluence
        URep = URep + (0.5 * repulsiveWeight * ((1/obstacle_distance) - (1/obstacle.distanceInfluence))^2);
    end
end
end

function [URepGrad] = potential_repulsiveSphereGrad(q)
global obstacles
global repulsiveWeight

URepGrad = zeros(2,1);
for i = 1:numel(obstacles)
    obstacle = obstacles(i);
    obstacle_distance = distance_toSphere(q, obstacle);
    obstacle_distanceGrad = distance_toSphereGrad(q, obstacle);
    if obstacle_distance <= obstacle.distanceInfluence
        URepGrad = URepGrad + repulsiveWeight * ((1/obstacle.distanceInfluence) - (1/obstacle_distance))*(1/obstacle_distance)^2 * obstacle_distanceGrad;
    end
end
end

%% Attractive Potential
function [UAttr] =  potential_attractive(q)
global goal;
global d_star;
global zeta;

dist = distance(q, goal);

if dist <= d_star
    UAttr = (1/2) * zeta * dist^2;
else
    UAttr = (dist * d_star * zeta) - (1/2) * zeta * d_star^2;
end
end

function [UAttrGrad] = potential_attractiveGrad(q)
global goal;
global d_star;
global zeta;

dist = distance(q, goal);
if dist <= d_star
    UAttrGrad = zeta * (q - goal);
else
    UAttrGrad = d_star * zeta * (q - goal) / dist ;
end
end

%% Total Potential
function [UGrad] = potential_totalGrad(q)
[URepGrad] = potential_repulsiveSphereGrad(q);
[UAttrGrad] = potential_attractiveGrad(q);

UGrad = UAttrGrad + URepGrad ;
end

%% Calculate Path
function [path] = potential_planner(start_)
global start;
global steps;
global epsilon;
global goal;

if isnan(start_)
    start_ = start;
end

path = zeros(2, steps);
path(:,1) = start_;

for step = 2:steps
    q = path(:,step - 1);
    
    UGrad = potential_totalGrad(path(:,step - 1));
    path(:,step) = path(:,step - 1) - epsilon * UGrad;
    dist_toGoal = distance(q, goal);
    if dist_toGoal < 0.005
        path = path(:,1:step);
        break;
    end
end
end

%% Helper Functions

%Calculate distance between 2 points
function [dist] = distance(q1, q2)
dist = norm(q1 - q2);
end

% Calculate distance between a point and a circle
function [dist] = distance_toSphere(q,sphere)
dist = distance(q, sphere.center) - sphere.radius;
end

% Calculate gradient of distance between point and circle
function [distGrad] = distance_toSphereGrad(q,sphere)
dist = distance_toSphere(q,sphere);
distGrad = (q - sphere.center)/dist;
end