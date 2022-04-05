function curr_fig = SingleSliceVisualization_ActivityMap(varargin)

% Visualization of the activity over slices around a fixed location
% The visualization shows three views: Axial, coronal, and sagittal views
% across the brain. A schematic insert helping to visualize the cut 
% lines are included.
%
% Input:
%     coord: (3,Ncoord) array, source space indicating the dipole locations
%                
%     Q: (1,ncoord) vector of non-negative entries, the dipole amplitudes 
%                
%     r0: (3,1) array, Cartesian coordinates of the cut planes
% 
% 
% Optional input: 
%   **Property**: 'Locator'. 
%   **Values**: 'haircross'(default),'circle','pinhead','none'
%     This property defines how the point r0 is marked in the plots
%
%   **Property**: 'Thickness'
%   **Values**: Real between 0 and 1. Default value 0.1
%     Specifies the thichness of the slice represented in the plot
%
%   **Property**: 'CoordSystem'
%   **Values**: 'MNI' (default) or 'Brainstorm'
%   *MNI*:        x (right), y (front), z (crown)
%   *Brainstorm*: x (front), y (left), z (crown)
% 
%
%  Output:
%     curr_fig: the figure handle
% 
% 
%  Usage:
%     SingleSliceVisualization(coord,Q,r0) 
%     SingleSliceVisualization(coord,Q,Y,r0,Property1,value1,Property2,value2,..)
%

%--------------------------------------------------------------------------
% CALLS TO: None
% Last updated: 7/22/2019
%--------------------------------------------------------------------------
%
n = length(varargin);
if n < 3
    disp('Not enough input arguments');
    return
end

% Necessary inputs

coord = varargin{1};
Q     = varargin{2};
r0    = varargin{3};

% Default values of the optional input parameters

RelThickness = 0.1;
locator = 'crosshair';
colorscale = 'single';
CoordSystem = 'MNI';

for j = 4:n
    if strcmp(varargin{j},'ColorScale')
        colorscale = varargin{j+1};
    end
    if strcmp(varargin{j},'Locator')
        locator = varargin{j+1};
    end
    if strcmp(varargin{j},'Thickness')
        RelThickness = varargin{j+1};
    end
    if strcmp(varargin{j},'CoordSystem')
        CoordSystem = varargin{j+1}; % 'MNI' (default) or 'Brainstorm'
    end
end

if strcmp(CoordSystem,'Brainstorm')
    % rotate the coordinate points in xy-plane by -pi/2
    Rotate = [0,-1,0;
              1, 0,0;
              0, 0,1];
    coord = Rotate*coord;
    r0    = Rotate*r0;
end

% Internal figure parameters:
% rel_width  : The width of the image realtive to the screen width
% rel_height : The height of the image relative to the screen height
% hor_posit  : Horizontal position of the lower left corner, relative to
%              lower left corner of the screen, units in screen width
% vert_posit : Vertical position of the lower left corner, relative to the
%              lower left corner of the screen, units in screen height
% locator    : String, defining how the selected position r0 in indicated.
%              Options: 'crosshair' (default) , 'circle','pinhead','none'
% colorscale : String, defines the color scale of the activity plot.
%              'single' (default): Four shades of blue + yellow 
%              'double': Five shades of blue + yellow, orange,red
% RelThickness = 0.1 (default), thichness of the slice 

rel_width  = 5/10;
rel_height = 6/10;
hor_posit  = 2/20;
vert_posit = 3/20;

% Generate the model for producing the explanatory insert figure, 
% consisting of a coarse triangulated geometric brain model. The (3,Nnodes) 
% matrix  Y contains the surface vertices, the (Nelem,3) matrix T is the 
% triangle topology. The model is based on a triangulated sphere that is
% read from a file.

load SphereTriangulation X T  % Generated with TriangulateSphere.m
[Y,T] = GetInsertModel(coord,X,T);

% Setting the window size, position, and background color
scrsz      = get(0,'ScreenSize');
width      = rel_width*scrsz(3);
height     = rel_height*scrsz(4);
curr_fig   = figure('Position',[hor_posit*scrsz(3),vert_posit*scrsz(4), width,height]);
background = [0,0,0]; % [0,0,0] = black background
set(gcf,'color',background);

Nnodes = size(Y,2); % Number of nodes in the surface mesh of the insert model
Nelem  = size(T,1);  % Number of elements in the surface mesh
Ncoord = size(coord,2); %Size of the source space
margin = 0.1; % Margin in the plots to avoid that figures get too crowded

xmin = min(Y(1,:));
xmax = max(Y(1,:));
ymin = min(Y(2,:));
ymax = max(Y(2,:));
zmin = min(Y(3,:));
zmax = max(Y(3,:));
 
% Inserts: Explain the slice locations 
% =====================================
 
 % Plotting a transparent head surface
 headcolor = [0.3,0.3,0.3];
% Position vectors for the geometric inserts
shift = 1/12;
pos_vec = [0,  0.3,0.6;
           3/4,3/4,3/4; 
           1/5,1/5,1/5;
           1/5,1/5,1/5];
 pos_vec(1,:) = pos_vec(1,:) + shift;

 % Axial view; use sagittal fig as an insert
pos_vecj = pos_vec(:,1)'; 
infig= axes('Position',pos_vecj);
 for j=1:size(T,1)
      Yj = [Y(:,T(j,1)),Y(:,T(j,2)),Y(:,T(j,3))];
      %Sagittal projection
     fill(Yj(2,:),Yj(3,:),headcolor,'EdgeColor',headcolor,'FaceAlpha',0.4)
     hold on
 end
axis('equal')
vv = axis;
plot([vv(1),vv(2)],[r0(3),r0(3)],'r-','LineWidth',3)
axis('off')
set(gca,'color',background);
hold off
text(vv(1),r0(3),'Back','FontSize',15,'Color','w','VerticalAlignment','bottom')

 % Coronal view; Use axial projection for explanation
 pos_vecj = pos_vec(:,2)'; 
infig= axes('Position',pos_vecj);
 for j=1:size(T,1)
      Yj = [Y(:,T(j,1)),Y(:,T(j,2)),Y(:,T(j,3))];
      %Sagittal projection
     fill(Yj(1,:),Yj(2,:),headcolor,'EdgeColor',headcolor,'FaceAlpha',0.4)
     hold on
 end
axis('equal')
vv = axis;
plot([vv(1),vv(2)],[r0(2),r0(2)],'r-','LineWidth',3)
axis('off')
set(gca,'color',background);
hold off
text(vv(1),r0(2),'Left','FontSize',15,'Color','w','VerticalAlignment','bottom')
 
% Sagittal view; Use axial projection for explanation
pos_vecj = pos_vec(:,3)'; 
infig= axes('Position',pos_vecj);
 for j=1:size(T,1)
      Yj = [Y(:,T(j,1)),Y(:,T(j,2)),Y(:,T(j,3))];
      %Sagittal projection
     fill(Yj(2,:),-Yj(1,:),headcolor,'EdgeColor',headcolor,'FaceAlpha',0.4)
     hold on
 end
axis('equal')
vv = axis;
add = 0.2*(vv(4)-vv(3));
plot([vv(1),vv(2)],[-r0(1),-r0(1)],'r-','LineWidth',3)
axis('off')
set(gca,'color',background);
hold off
text(vv(1),-r0(1),'Back','FontSize',15,'Color','w','VerticalAlignment','bottom')
  
 % Activity plots
 % ==============
 
if strcmp(colorscale,'single')
    % Classifying the dipoles in five strength categories
    % Here we use a linear color scale
    Qmax = max(abs(Q));
    ncolors = 5;
    Colors = cell(5,3);
    Colors{1,1} = find(Q<0.3*Qmax);
    Colors{2,1} = find(Q>=0.3*Qmax & Q <0.5*Qmax);
    Colors{3,1} = find(Q>=0.5*Qmax & Q <0.7*Qmax);
    Colors{4,1} = find(Q>=0.7*Qmax & Q <0.9*Qmax);
    Colors{5,1} = find(Q>=0.9*Qmax);
    % Blue scale
    Colors{1,2} = 1/255*[0,0,145];   % Dark Blue
    Colors{2,2} = 1/255*[0,0,190];   % Medium Blue
    Colors{3,2} = 1/255*[0,0,255];   % Blue
    Colors{4,2}= 1/255*[0,191,255];  % Deep Sky Blue
    Colors{5,2} = 1/255*[255,255,0]; % Yellow
    % Marker sizes
    Colors{1,3} = 5;
    Colors{2,3} = 10;
    Colors{3,3} = 10;
    Colors{4,3} = 10;
    Colors{5,3} = 10;
elseif strcmp(colorscale,'double')
    
    
    
    
    
    
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
rgb_grey = [128, 128, 128];
Colors{1,2} = 1/255*rgb_grey; 
Colors{2,2} = 1/255*rgb_grey; 
Colors{3,2} = 1/255*rgb_grey; 
Colors{4,2} = 1/255*rgb_grey; 
% Colors{4,2} = 1/255*[50,50,50]; 
% Colors{5,2} = 1/255*[50,50,50]; 
% Colors{6,2} = 1/255*[50,50,50]; 
% Colors{4,2} = 1/255*[0,75,125]; 
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
end
background = 1/255*[0,0,0];
%     % Classifying the dipoles in five strength categories
%     % The color scale is piecewise linear in logarithmic scale
%     ncolors = 8; 
%     eta = 0.90;
%     mu = 0.35;
%     Qmin = min(Q);
%     Qmax = max(Q);
%     Qsort = sort(Q,'ascend');
%     icut = round(eta*length(Q));
%     % Moving min(Q) up so that the low outliers do not distort the color
%     % scale
%     ilow = round(mu*length(Q));
%     Qmin = Qsort(ilow);
%     Qcut = Qsort(icut);
%     logQ = log10(Q);
%     logQlev1 = linspace(log10(Qmin),log10(Qcut),5);
%     logQlev2 = linspace(log10(Qcut),log10(Qmax),5);
%     logQlev   = [logQlev1, logQlev2(2:end)];
%     Colors = cell(8,3);
%     Colors{1,1} = find(logQ<logQlev(2));
%     Colors{2,1} = find(logQ>=logQlev(2) & logQ <logQlev(3));
%     Colors{3,1} = find(logQ>=logQlev(3) & logQ <logQlev(4));
%     Colors{4,1} = find(logQ>=logQlev(4) & logQ <logQlev(5));
%     Colors{5,1} = find(logQ>=logQlev(5) & logQ <logQlev(6));
%     Colors{6,1} = find(logQ>=logQlev(6) & logQ <logQlev(7));
%     Colors{7,1} = find(logQ>=logQlev(7) & logQ <logQlev(8));
%     Colors{8,1} = find(logQ>=logQlev(8));
%     % Double color scale
%     Colors{1,2} = 1/255*[0,0,145];   % Dark Blue
%     Colors{2,2} = 1/255*[0,0,190];   % Medium Blue
%     Colors{3,2} = 1/255*[0,0,255];   % Blue
%     Colors{4,2}= 1/255*[0,125,255];  % Deep Sky Blue
%     Colors{5,2} = 1/255*[0,255,255]; % Cyan
%     Colors{6,2} = 1/255*[255,255,0]; % Yellow
%     Colors{7,2} = 1/255*[255,165,0]; % Orange
%     Colors{8,2} = 1/255*[255,10,10]; % Red 
%     % Marker sizes
%     Colors{1,3} = 10;
%     Colors{2,3} = 10;
%     Colors{3,3} = 10;
%     Colors{4,3} = 10;
%     Colors{5,3} = 10;
%     Colors{6,3} = 10;
%     Colors{7,3} = 10;
%     Colors{8,3} = 10;
% end
 
pos_vec = [0,   0.3, 0.6;
           0,   0,   0; 
           0.3, 0.3, 0.4;
           7/10,7/10,7/10];

zmin = min(coord(3,:));
zmax = max(coord(3,:));
xmin = min(coord(1,:));
xmax = max(coord(1,:));
ymin = min(coord(2,:));
ymax = max(coord(2,:));
 
 % Axial view
 
Thickness= RelThickness*(zmax - zmin);
slice = find(coord(3,:)>=r0(3) - Thickness/2 & coord(3,:)<=r0(3) + Thickness/2);
pos_vecj = pos_vec(:,1)'; 
actfig= axes('Position',pos_vecj);
 
 for k = 1:ncolors
        Ij = intersect(slice,Colors{k,1});
        plot(coord(1,Ij),coord(2,Ij),'.','color',Colors{k,2},'MarkerSize',Colors{k,3});
        hold on
 end
 dx = xmax - xmin;
 dy = ymax - ymin;
 xlow = xmin - margin*dx;
 xhigh = xmax + margin*dx;
 ylow = ymin - margin*dy;
 yhigh = ymax + margin*dy;
 axis([xlow,xhigh,ylow,yhigh]);
 axis('equal')
 axis('off')
 
 % Marking the true sources

nsource = size(r0,2);


if strcmp(locator,'pinhead')
    % Use pinhead markers for true sources
    for j = 1:nsource
        zsource = r0(3,j);
        if strcmp(colorscale,'single')
           plot(actfig,[r0(1,j),r0(1,j)],[r0(2,j),ymax +  0.03*(ymax-ymin)],'-r','LineWidth',1)
           plot(actfig,r0(1,j),ymax +  0.03*(ymax-ymin),'rs','MarkerSize',8,'MarkerFaceColor','r')
        elseif strcmp(colorscale,'double')
           plot(actfig,[r0(1,j),r0(1,j)],[r0(2,j),ymax +  0.03*(ymax-ymin)],'-g','LineWidth',1)
           plot(actfig,r0(1,j),ymax +  0.03*(ymax-ymin),'gs','MarkerSize',8,'MarkerFaceColor','g')
        end
    end
elseif strcmp(locator,'crosshair')
    % Use crosshair markers
    tx = 0.03*(xmax-xmin);
    ty = 0.03*(ymax-ymin);
    for j = 1:nsource
        zsource = r0(3,j);
        plot(actfig,[xmin-tx,xmax+tx],[r0(2,j),r0(2,j)],'r-','LineWidth',1)
        plot(actfig,[r0(1,j),r0(1,j)],[ymin-ty,ymax+ty],'r-','LineWidth',1)
    end
 elseif strcmp(locator,'circle')
    th = linspace(0,2*pi,20);
    rho = 0.12*(xmax-xmin);
    for j = 1:nsource
        plot(actfig,r0(1,j) +rho*cos(th),r0(2,j) + rho*sin(th),'r-','LineWidth',3)
    end
end
hold off

 % Coronal view
 
 Thickness= RelThickness*(ymax - ymin);
 slice = find(coord(2,:)>=r0(2) - Thickness/2 & coord(2,:)<=r0(2) + Thickness/2);
  pos_vecj = pos_vec(:,2)'; 
 actfig= axes('Position',pos_vecj);
 
 for k = 1:ncolors
        Ij = intersect(slice,Colors{k,1});
        plot(coord(1,Ij),coord(3,Ij),'.','color',Colors{k,2},'MarkerSize',Colors{k,3});
        hold on
 end

 dz = zmax - zmin;
 zlow = zmin - margin*dz;
 zhigh = zmax + margin*dz;
 axis([xlow,xhigh,zlow,zhigh]);
 axis('equal')
 axis('off')
 
 
nsource = size(r0,2);

if strcmp(locator,'pinhead')
    % Use pinhead markers for true sources
    for j = 1:nsource
        ysource = r0(2,j);
        if strcmp(colorscale,'single')
           plot(actfig,[r0(1,j),r0(1,j)],[r0(3,j),zmax +  0.03*(zmax-zmin)],'r-','LineWidth',1)
           plot(actfig,r0(1,j),zmax +  0.03*(zmax-zmin),'rs','MarkerSize',8,'MarkerFaceColor','r')
        elseif strcmp(colorscale,'double')
           plot(actfig,[r0(1,j),r0(1,j)],[r0(3,j),zmax +  0.03*(zmax-zmin)],'g-','LineWidth',1)
           plot(actfig,r0(1,j),zmax +  0.03*(zmax-zmin),'gs','MarkerSize',8,'MarkerFaceColor','g')
        end
    end
    
elseif strcmp(locator,'crosshair')
    tx = 0.03*(xmax-xmin);
    tz = 0.03*(zmax-zmin);
    % Use crosshair markers
    for j = 1:nsource
        ysource = r0(2,j);
        plot(actfig,[xmin-tx,xmax+tx],[r0(3,j),r0(3,j)],'r-','LineWidth',1)
        plot(actfig,[r0(1,j),r0(1,j)],[zmin-tz,zmax+tz],'r-','LineWidth',1)
    end
elseif strcmp(locator,'circle')
    th = linspace(0,2*pi,20);
    rho = 0.12*(xmax-xmin);
    for j = 1:nsource
        plot(actfig,r0(1,j) +rho*cos(th),r0(3,j) + rho*sin(th),'r-','LineWidth',3)
    end
end
hold off

 % Sagittal view
 
  Thickness= RelThickness*(xmax - xmin);
 slice = find(coord(1,:)>=r0(1) - Thickness/2 & coord(1,:)<=r0(1) + Thickness/2);
  pos_vecj = pos_vec(:,3)'; 
 actfig= axes('Position',pos_vecj);
 
 for k = 1:ncolors
        Ij = intersect(slice,Colors{k,1});
       plot(coord(2,Ij),coord(3,Ij),'.','color',Colors{k,2},'MarkerSize',Colors{k,3});
        hold on
 end
 axis([ylow,yhigh,zlow,zhigh]);
 axis('equal')
 axis('off')
 nsource = size(r0,2);

if strcmp(locator,'pinhead')
    % Use pinhead markers for true sources
    for j = 1:nsource
        xsource = r0(1,j);
        if strcmp(colorscale,'single')
           plot(actfig,[r0(2,j),r0(2,j)],[r0(3,j),zmax + 0.03*(zmax-xmin)],'r-','LineWidth',1)
           plot(actfig,r0(2,j),zmax + 0.03*(zmax-xmin),'rs','MarkerSize',8,'MarkerFaceColor','r')
        elseif strcmp(colorscale,'double')
           plot(actfig,[r0(2,j),r0(2,j)],[r0(3,j),zmax + 0.03*(zmax-xmin)],'g-','LineWidth',1)
           plot(actfig,r0(2,j),zmax + 0.03*(zmax-xmin),'gs','MarkerSize',8,'MarkerFaceColor','g')
        end
        
    end
    
elseif strcmp(locator,'crosshair')
    % Use crosshair markers
    ty = 0.03*(ymax-ymin);
    tz = 0.03*(zmax-zmin);
    for j = 1:nsource
        xsource = r0(1,j);
        plot(actfig,[ymin-ty,ymax+ty],[r0(3,j),r0(3,j)],'r-','LineWidth',1)
        plot(actfig,[r0(2,j),r0(2,j)],[zmin-tz,zmax+tz],'r-','LineWidth',1)
    end
 elseif strcmp(locator,'circle')
    th = linspace(0,2*pi,20);
    rho = 0.1*(ymax-ymin);
    for j = 1:nsource
         plot(actfig,r0(2,j) +rho*sin(th),r0(3,j) + rho*cos(th),'r-','LineWidth',3)
    end
end

hold off

% 3. Explanatory coordinate frame

pos_vec = [0,   1/3, 2/3;
           6/10,6/10,6/10; 
           1/20,1/20,1/20;
           1/10,1/10,1/10];
 shift = 0.01;
 pos_vec(1,:) = pos_vec(1,:) + shift;
 
 pos_vecj = pos_vec(:,1)'; 
 actfig= axes('Position',pos_vecj);
 plot([0,0],[0,1],'w-','Linewidth',2)
 hold on
 plot([0,1],[0,0],'w-','Linewidth',2)
 set(gca,'color',background);
 text(0.1,1,'Front','Color','w','HorizontalAlignment','left','VerticalAlignment','cap','FontSize',15)
 text(1,0,'Right','Color','w','HorizontalAlignment','right','VerticalAlignment','top','FontSize',15)
 hold off
  
 pos_vecj = pos_vec(:,2)'; 
 actfig= axes('Position',pos_vecj);
 plot([0,0],[0,1],'w-','Linewidth',2)
 hold on
 plot([0,1],[0,0],'w-','Linewidth',2)
 set(gca,'color',background);
 text(0.1,1,'Crown','Color','w','HorizontalAlignment','left','VerticalAlignment','cap','FontSize',15)
 text(1,0,'Right','Color','w','HorizontalAlignment','right','VerticalAlignment','top','FontSize',15)
 hold off
  
 pos_vecj = pos_vec(:,3)'; 
 actfig= axes('Position',pos_vecj);
 plot([0,0],[0,1],'w-','Linewidth',2)
 hold on
 plot([0,1],[0,0],'w-','Linewidth',2)
 set(gca,'color',background);
 text(0.1,1,'Crown','Color','w','HorizontalAlignment','left','VerticalAlignment','cap','FontSize',15)
 text(1,0,'Front','Color','w','HorizontalAlignment','right','VerticalAlignment','top','FontSize',15)
 hold off
 
 % Turning of the InvertHardcopy to allow printing of the figure with the
 % current background color

curr_fig.InvertHardcopy ='off';
 
  
