function [APChol, has_normal_vec] = BuildAnatomicalPrior(coord, normals, tiny)
% This function generates the anatomical prior that favors dipoles in the preferred direction. 
% 
% Input:
%      coord: (3,N) array, coordinates of the dipoles in the grid
% 
%      normals: (3,N) array, normal vectors of the dipoles in the grid
% 
%      tiny: scalar, relative variance of the components of the dipoles orthogonal to the preferred direction (default value: 0.05)
%
% Output:
%       APChol: (3*N,3*N) array, the Cholesky factor of the anatomical prior covariance 
%
% Usage: 
%       APChol = BuildAnatomicalPrior(coord,normals,tiny)
%       BuildAnatomicalPrior(coord,normals) is equivalent to BuildAnatomicalPrior(coord,normals,0.05)
% 


% Setting default values
if nargin == 2, tiny = 0.05; end

% Reading the total number N of grid points in the brain model
N = size(coord, 2);   

% Costructing a local orthonormal frame at each grid point
U1 = NaN(3, N);
U2 = NaN(3, N);
U3 = NaN(3, N);
has_normal_vec = NaN(N, 1);
E = eye(3);
for ell = 1:N
    u3 = normals(:,ell);
    if sum(u3) == 0
      has_normal_vec(ell) = false;
    else
      has_normal_vec(ell) = true;
    end
    [u3sort,Isort] = sort(abs(u3),'ascend');   % q3(3) largest
    u3sort = u3sort.*sign(u3(Isort));
    Pell = E(Isort,:);
    u1 = [u3sort(3);0;-u3sort(1)];
    u1 = Pell'*u1;
    u2 = [u3(2)*u1(3) - u1(2)*u3(3);u3(3)*u1(1)-u1(3)*u3(1);u3(1)*u1(2)-u1(1)*u3(2)];
    U1(:,ell) = (1/norm(u1))*u1;
    U2(:,ell) = (1/norm(u2))*u2;
    U3(:,ell) = (1/norm(u3))*u3;
end

% Computing the local anatomical prior covariance matrices   
C = sparse(3*N,3*N);
for j = 1:N
  if has_normal_vec(j) 
    C(((j-1)*3+1):j*3,((j-1)*3+1):3*j) = tiny*(U1(:,j)*U1(:,j)' + U2(:,j)*U2(:,j)') ...
            + U3(:,j)*U3(:,j)';
  else
    C(((j-1)*3+1):j*3,((j-1)*3+1):3*j) = (1+2*tiny)/3 * eye(3);
  end
end

APChol = chol(C);
