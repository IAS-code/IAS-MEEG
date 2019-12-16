function [Y,T] = GetInsertModel(coord,X,T)


% This program generates a simple surface model enclosing the source space
% defined by the coordinates  given by the input matrix coord.
%
% Input: coord   - (3,N) matrix, defining the source space
%        X       - (3,n) matrix, nodes on the sphere's surface
%        T       - (n_elem,3) topology matrix
%
% Output: Y      - (3,ny) matrix, coordinates of the nodes
%         T      - (n_elem,3) matrix, topology defining the triangles
%--------------------------------------------------------------------------
% CALLS TO: None
% 2/20/2015 - Roma. Modified 11/26/19
%--------------------------------------------------------------------------

N = size(coord,2);        % Size of the source space

% Geometric center point of the source space

x_min = min(coord(1,:));
x_max = max(coord(1,:));
y_min = min(coord(2,:));
y_max = max(coord(2,:));
z_min = min(coord(3,:));
z_max = max(coord(3,:));
c = 0.5*[x_min+x_max;y_min+y_max;z_min+z_max];

% Computing the surface triangularization of the brain

coord = coord - c*ones(1,N); % Center the coordinate system

n = size(X,2);
Y = NaN(3,n);
% To compute the projection, select source points that
% are in a cone around the direction, cos(angle)>coscone
coscone   = 0.95;
normcoord = sqrt(sum(coord.^2,1));
inflate   = 1.02;
for j = 1:n
    alpha = X(:,j); 
    % Restricting the coordinates involved into a cone
    p = alpha'*coord;
    cc = p./normcoord;
    I = find(cc>coscone);
    if length(I) == 0
        pj = 0;
    else
       [pj,ij] = max(p(I));
    end
    Y(:,j) = inflate*pj*alpha;
end

% Shifting the points around a new center

Y = c*ones(1,size(Y,2)) + Y;

%figure(1)
%for j=1:size(T,1)
%    Yj = [Y(:,T(j,1)),Y(:,T(j,2)),Y(:,T(j,3))];
%    fill3(Yj(1,:),Yj(2,:),Yj(3,:),[0.4,0.4,0.8],'FaceAlpha',0.7);
%    axis('equal')
%    hold on
%end
%hold off

 
 
