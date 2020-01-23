function [id_a, id_c, id_s] = SlicedVisualization_ActivityMap(coord, Q)
%
% This is a development program for visualizing the activity map
%
% Input:
%      coord: (3,N) array, coordinates of the dipoles in the grid
%
%      Q: (3*N, 1) vector, reconstructed dipole moments at a specified time point
%      
% Usage:
%      SlicedVisualization_ActivityMap_Gif(coord, Q_est);


%% Classifying the dipoles in five categories according to their strength 
% The color scale is piecewise linear in logarithmic scale
ncolors = 8; 
eta = 0.9;
mu = 0.35;
Qmax = max(Q);
Qsort = sort(Q,'ascend');
icut = round(eta*length(Q));
% Moving min(Q) up so that the low outliers do not distort the color
% scale
ilow = round(mu*length(Q));
Qmin = Qsort(ilow);
Qcut = Qsort(icut);
logQ = log10(Q);
logQlev1 = linspace(log10(Qmin),log10(Qcut),5);
logQlev2 = linspace(log10(Qcut),log10(Qmax),5);
logQlev   = [logQlev1, logQlev2(2:end)];
Colors = cell(8,3);
Colors{1,1} = find(logQ<logQlev(2));
Colors{2,1} = find(logQ>=logQlev(2) & logQ <logQlev(3));
Colors{3,1} = find(logQ>=logQlev(3) & logQ <logQlev(4));
Colors{4,1} = find(logQ>=logQlev(4) & logQ <logQlev(5));
Colors{5,1} = find(logQ>=logQlev(5) & logQ <logQlev(6));
Colors{6,1} = find(logQ>=logQlev(6) & logQ <logQlev(7));
Colors{7,1} = find(logQ>=logQlev(7) & logQ <logQlev(8));
Colors{8,1} = find(logQ>=logQlev(8));
% Double color scale
Colors{1,2} = 1/255*[50,50,50]; 
Colors{2,2} = 1/255*[50,50,50]; 
Colors{3,2} = 1/255*[50,50,50]; 
% Colors{4,2} = 1/255*[50,50,50]; 
% Colors{5,2} = 1/255*[50,50,50]; 
% Colors{6,2} = 1/255*[50,50,50]; 
Colors{4,2} = 1/255*[0,75,125]; 
Colors{5,2} = 1/255*[0,125,125]; 
Colors{6,2} = 1/255*[255,255,0]; 
Colors{7,2} = 1/255*[255,165,0];
Colors{8,2} = 1/255*[255,10,10];
% Marker sizes
Colors{1,3} = 10;
Colors{2,3} = 10;
Colors{3,3} = 10;
Colors{4,3} = 10;
Colors{5,3} = 10;
Colors{6,3} = 10;
Colors{7,3} = 10;
Colors{8,3} = 10;
background = 1/255*[0,0,0];

%% Setting slices
zmin = min(coord(3,:));
zmax = max(coord(3,:));
xmin = min(coord(1,:));
xmax = max(coord(1,:));
ymin = min(coord(2,:));
ymax = max(coord(2,:));
nslice = 10;
Slices = cell(nslice,3); % Columns: axial, coronal,sagittal views
                         % Rows:    slice by slice

% Axial slicing: From bottom to top
z = linspace(zmin,zmax,nslice+1);
count = 0;
for j = 1:nslice
    count = count + 1;
    Slices{count} = find(coord(3,:)>=z(j) & coord(3,:)<=z(j+1));
end

% Sagittal slicing: From left to right
y = linspace(ymin,ymax,nslice+1);
for j = 1:nslice
    count = count + 1;
    Slices{count}= find(coord(2,:)>=y(j) & coord(2,:)<=y(j+1));
end

% Coronal slicing: From back to front
x = linspace(xmin,xmax,nslice+1);
for j = 1:nslice
    count = count + 1;
    Slices{count} = find(coord(1,:)>=x(j) & coord(1,:)<=x(j+1));
end

%% Plotting
v1 = (0:round(nslice/2)-1);
v1 =2/nslice* [v1,v1];
v2 = [0.5*ones(1,round(nslice/2)),zeros(1,round(nslice/2))];

% Axial slices
scrsz = get(0,'ScreenSize');
id_a = figure('Name','Axial slices','Position',[scrsz(3)/20,scrsz(4)/20, 18*scrsz(3)/20, 18*scrsz(4)/20]);
set(id_a,'color',background);
for j = 1:nslice
    pos_vec1 = [v1(j),v2(j),2/nslice,1/2];
    AxialFig1=axes('Position',pos_vec1);
    for k = 1:ncolors
        Ij = intersect(Slices{j,1},Colors{k,1});
        plot(AxialFig1,coord(1,Ij),coord(2,Ij),'.','color',Colors{k,2},'MarkerSize',Colors{k,3});
        hold on
    end
    axis(AxialFig1,[xmin,xmax,ymin,ymax]);
    axis(AxialFig1,'equal')
    axis(AxialFig1,'off')
    % Plot axial sections
    pos_vec2 = [v1(j)+0.15,v2(j),2/nslice*0.2,1/2*0.2];
    AxialFig2=axes('Position',pos_vec2);
    plot(AxialFig2,coord(1,:),coord(3,:),'.','color',Colors{1,2},'MarkerSize',Colors{1,3});
    hold on
    Ij = Slices{j,1};
    plot(AxialFig2,coord(1,Ij),coord(3,Ij),'.','color',Colors{8,2},'MarkerSize',Colors{8,3});
    axis(AxialFig2,[xmin,xmax,zmin,zmax]);
    axis(AxialFig2,'equal')
    axis(AxialFig2,'off')
end
hold off


% Sagittal slices
id_c = figure('Name','Sagittal slices','Position',[scrsz(3)/20,scrsz(4)/20, 18*scrsz(3)/20, 18*scrsz(4)/20]);
set(id_c,'color',background);
for j = 1:nslice
    pos_vec1 = [v1(j),v2(j),2/nslice,1/2];
    SagittalFig1=axes('Position',pos_vec1);
    for k = 1:ncolors
        Ij = intersect(Slices{j,2},Colors{k,1});
        plot(SagittalFig1,coord(1,Ij),coord(3,Ij),'.','color',Colors{k,2},'MarkerSize',Colors{k,3});
        hold on
    end
    axis(SagittalFig1,[xmin,xmax,zmin,zmax]);
    axis(SagittalFig1,'equal')
    axis(SagittalFig1,'off')
    % Plot Sagittal sections
    pos_vec2 = [v1(j)+0.15,v2(j),2/nslice*0.2,1/2*0.2];
    SagittalFig2=axes('Position',pos_vec2);
    plot(SagittalFig2,coord(1,:),coord(2,:),'.','color',Colors{1,2},'MarkerSize',Colors{1,3});
    hold on
    Ij = Slices{j,2};
    plot(SagittalFig2,coord(1,Ij),coord(2,Ij),'.','color',Colors{8,2},'MarkerSize',Colors{1,3});
    axis(SagittalFig2,[xmin,xmax,ymin,ymax]);
    axis(SagittalFig2,'equal')
    axis(SagittalFig2,'off')
end
hold off

% Coronal slices
id_s = figure('Name','Coronal slices','Position',[scrsz(3)/20,scrsz(4)/20, 18*scrsz(3)/20, 18*scrsz(4)/20]);
set(id_s,'color',background);
for j = 1:nslice
    pos_vec1 = [v1(j),v2(j),2/nslice,1/2];
    CoronalFig1=axes('Position',pos_vec1);
    for k = 1:ncolors
        Ij = intersect(Slices{j,3},Colors{k,1});
        plot(CoronalFig1,coord(3,Ij),coord(2,Ij),'.','color',Colors{k,2},'MarkerSize',Colors{k,3});
        hold on
    end
    axis(CoronalFig1,[zmin,zmax,ymin,ymax]);
    axis(CoronalFig1,'equal')
    axis(CoronalFig1,'off')
    % Plot Coronal sections
    pos_vec2 = [v1(j)+0.15,v2(j),2/nslice*0.2,1/2*0.2];
    CoronalFig2=axes('Position',pos_vec2);
    plot(CoronalFig2,coord(1,:),coord(3,:),'.','color',Colors{1,2},'MarkerSize',Colors{1,3});
    hold on
    Ij = Slices{j,3};
    plot(CoronalFig2,coord(1,Ij),coord(3,Ij),'.','color',Colors{8,2},'MarkerSize',Colors{1,3});
    axis(CoronalFig2,[xmin,xmax,zmin,zmax]);
    axis(CoronalFig2,'equal')
    axis(CoronalFig2,'off')
end
hold off

