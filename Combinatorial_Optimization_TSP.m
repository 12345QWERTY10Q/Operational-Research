%%
% Bivariate Integer Programming
% Combinatorial Optimization TSP
%%
% 问题表示
% 1.生成所有可能的行程，意味着经过所有不同停留点对组。
% 2.计算每次行程的距离。
% 3.要最小化的代价函数是行程中每次旅行的旅行距离之和。
% 4.决策变量是与每个行程相关联的二元变量，其中每个 1 表示环程中存在的一次行程，每个 0 表示不在环程中的一次行程。
% 5.为确保环程包括每个经停留点，应加入这样一个线性约束：每个停留点都正好涉及两段行程。这意味着一段行程到达该停留点，一段行程离开该停留点。
clc;
clear all;
%%
% Generating stop
load('usborder.mat','x','y','xx','yy');
rng(3,'twister') % Makes stops in Maine & Florida, and is reproducible
nStops = 200; % You can use any number, but the problem size scales as N^2
stopsLon = zeros(nStops,1); % Allocate x-coordinates of nStops
stopsLat = stopsLon; % Allocate y-coordinates
n = 1;
while (n <= nStops)
    xp = rand*1.5;
    yp = rand;
    if inpolygon(xp,yp,x,y) % Test if inside the border
        stopsLon(n) = xp;
        stopsLat(n) = yp;
        n = n+1;
    end
end
%%
% Calculating the distance between stops
idxs = nchoosek(1:nStops,2);
dist = hypot(stopsLat(idxs(:,1)) - stopsLat(idxs(:,2)), ...
             stopsLon(idxs(:,1)) - stopsLon(idxs(:,2)));
lendist = length(dist);

%%
% Creating diagrams and draw maps
G = graph(idxs(:,1),idxs(:,2));
figure
hGraph = plot(G,'XData',stopsLon,'YData',stopsLat,'LineStyle','none','NodeLabel',{});
hold on
% Drawing the outside border
plot(x,y,'r-')
hold off
%%
% Creating variables and problems
tsp = optimproblem;
trips = optimvar('trips',lendist,1,'Type','integer','LowerBound',0,'UpperBound',1);
tsp.Objective = dist'*trips;
%%
% Constraints
constr2trips = optimconstr(nStops,1);
for stop = 1:nStops
    whichIdxs = outedges(G,stop); % Identify trips associated with the stop
    constr2trips(stop) = sum(trips(whichIdxs)) == 2;
end
tsp.Constraints.constr2trips = constr2trips;
%%
% Solving the initial problem
opts = optimoptions('intlinprog','Display','off');
tspsol = solve(tsp,'options',opts);
%%
% Visualizing resolution
tspsol.trips = logical(round(tspsol.trips));
Gsol = graph(idxs(tspsol.trips,1),idxs(tspsol.trips,2),[],numnodes(G));
hold on
highlight(hGraph,Gsol,'LineStyle','-')
title('Solution with Subtours')
%%
% Subloop path constraint
tourIdxs = conncomp(Gsol);
numtours = max(tourIdxs); % Number of subtours
fprintf('# of subtours: %d\n',numtours);
% Index of added constraints for subtours
k = 1;
while numtours > 1 % Repeat until there is just one subtour
    % Add the subtour constraints
    for ii = 1:numtours
        inSubTour = (tourIdxs == ii); % Edges in current subtour
        a = all(inSubTour(idxs),2); % Complete graph indices with both ends in subtour
        constrname = "subtourconstr" + num2str(k);
        tsp.Constraints.(constrname) = sum(trips(a)) <= (nnz(inSubTour) - 1);
        k = k + 1;        
    end
    
    % Try to optimize again
    [tspsol,fval,exitflag,output] = solve(tsp,'options',opts);
    tspsol.trips = logical(round(tspsol.trips));
    Gsol = graph(idxs(tspsol.trips,1),idxs(tspsol.trips,2),[],numnodes(G));
    % Gsol = graph(idxs(tspsol.trips,1),idxs(tspsol.trips,2)); % Also works in most cases
    
    % Plot new solution
    hGraph.LineStyle = 'none'; % Remove the previous highlighted path
    highlight(hGraph,Gsol,'LineStyle','-')
    drawnow

    % How many subtours this time?
    tourIdxs = conncomp(Gsol);
    numtours = max(tourIdxs); % Number of subtours
    fprintf('# of subtours: %d\n',numtours)    
end
title('Solution with Subtours Eliminated');
hold off
%%
% Checking the quality of solution
fprintf('Absolute Clearance:')
disp(output.absolutegap)
