function [Q,diagnostics] = IAS_algorithm(LF,LF_scaling,APChol,B,B_scaling,sigma,theta_star,eta,n_outer,maxit,theta_tol)
%
% This function solves the MEG inverse problem by the Iterative Alternating Sequential (IAS) algorithm. 
% The IAS algorithm is based on an iterative scheme that alternatively updates the dipole
% moments Q by solving a linear least squares problem using a priorconditioned CGLS algorithm with sutable 
% stopping condition and updating the hyperparameter theta by an explicit formula.
%
% Input:
%     LF: (M,3*N) array, the lead field matrix where M is the number of channels and N is the number of dipoles
% 
%     LF_scaling: scalar to be used to scale the lead field (output of the function SetParameters)
% 
%     APChol: (3*N,3*N) array, the Cholesky factor of the anatomical prior covariance matrix
% 
%     B: (M,T) array, a set of data of length T
% 
%     B_scaling: scalar, to be used for scaling the data (output of the function SetParameters)
% 
%     sigma: scalar, standard deviation of the scaled noise (output of the function SetParameters)
% 
%     theta_star: N-vector, computed theta_star values (output of the function SetParameters)
% 
%     eta: scalar, parameter for selecting the focality of the reconstructed activity (default value: 0.01; choose 0.001 for focal sources)
% 
%     n_outer: scalar, maximum number of iterations of the outer loop (default value: 30)
% 
%     maxit: scalar, maximum number of iterations of the inner loop (default value: 120)
% 
%     theta_tol: scalar, limit of the relative change of theta in the outer loop (default value: 0.001)
%
% 
% Output:
%      Q: (3*N,T) array, reconstructed dipole moments
% 
%      diagnostics: [n_outer+1,T] matrix containing the number of inner iterations of
%      each outer iteration step. Zeros in the matrix indicate that
%      convergence was reached. The last row returns the relative changes in theta at the last outer iteration. 
% 
% Usage:
%     [Q,diagnostics] = IAS_algorithm(LF,LF_scaling,APChol,B,B_scaling,sigma,theta_star,eta,n_outer,maxit,theta_tol)
% 
%     IAS_algorithm(LF,LF_scaling,APChol,B,B_scaling,sigma,theta_star) is equivalent to IAS_algorithm(LF,LF_scaling,APChol,B,B_scaling,sigma,theta_star,0.01,30,120,0.001)
% 
%     IAS_algorithm(LF,LF_scaling,APChol,B,B_scaling,sigma,theta_star,eta) is equivalent to IAS_algorithm(LF,LF_scaling,APChol,B,B_scaling,sigma,theta_star,eta,30,120,0.001)
% 
%     IAS_algorithm(LF,LF_scaling,APChol,B,B_scaling,sigma,theta_star,eta,n_outer,maxit) is equivalent to IAS_algorithm(LF,LF_scaling,APChol,B,B_scaling,sigma,theta_star,eta,n_outer,maxit,0.001)
%
%---------------------------------------------------------------
% References:
%
% D. Calvetti, A. Pascarella, F. Pitolli, E. Somersalo, B. Vantaggi
% A hierarchical Krylov-Bayes iterative inverse solver for MEG with physiological preconditioning
% Inverse Problems, 31 (12) (2015) 12500
% 
% D. Calvetti, A. Pascarella, F. Pitolli, E. Somersalo, B. Vantaggi
% Brain Activity Mapping from MEG Data via a Hierarchical Bayesian Algorithm with Automatic Depth Weighting,
% Brain topography, 1-31, (2018) 
%----------------------------------------------------------------
% Version: July, 2017
%----------------------------------------------------------------

%% Input

% Setting default values
if nargin == 7, eta = 0.01; n_outer = 30; maxit = 120; theta_tol = 0.001; end
if nargin == 8,             n_outer = 30; maxit = 120; theta_tol = 0.001; end
if nargin == 9,                           maxit = 120; theta_tol = 0.001; end
if nargin == 10,                                       theta_tol = 0.001; end

% Reading the number of channels M and the number of dipoles N
[M,N3] = size(LF);
N = N3/3;

% Target discrepancy of the whitened problem 
discr = sqrt(M); 

% Scaling to whiten the noise
A  = (1/sigma)*LF_scaling*LF;
BB = (1/sigma)*B_scaling*B;
          
% Whitening the unknown dipole moments with anatomical prior
A = A*APChol';

%% Time loop

% Setting the time interval of the protocol
nT = size(B,2);

diagnostics = zeros(n_outer+1,nT);

% Inizializing theta and Q
theta = theta_star;
Q = NaN(3*N,nT);

% Starting the time loop 
disp('IAS algorithm: Running time loop ')
for jj = 1:nT
     disp(['time = ', num2str(jj), ' (max. time = ',num2str(nT),')']);
	 % Extracting data at time jj
	 b = BB(:,jj);

	 % Solving the inverse problem at time jj using IAS algorithm
	 w = zeros(3*N,1);	 
	 iterate = 'yes';
	 outer_iteration='yes';
	 outer_count=0;
	 
	 % Outer loop for updating the hyperparameter theta
	 while (outer_count<n_outer)&strcmp(outer_iteration,'yes')
		 outer_count=outer_count+1;
		 Aj = A*spdiags(kron(sqrt(theta),ones(1,3))',0,3*N,3*N);       
		 % Step 1: Updating the dipole moments q using PCGLS
		 d = b - Aj*w;
		 r = Aj'*d;
		 p = r;
		 aux0 = r'*r;
		 j = 0;
		 iterate = 'yes';    

		% Inner loop for updating w
		 while (j<=maxit)&strcmp(iterate,'yes')
			 j = j + 1;
			 y = Aj*p;
			 alpha = aux0/(y'*y);
			 w = w + alpha*p;       % new iterate
			 d = d - alpha*y;       % new discrepancy
			 % Checking the discrepancy condition
			 if norm(d) < discr
			    % The discrepancy is below the threshold; stop the inner iteration
				iterate = 'no';
			 end
			 r = Aj'*d;
			 aux1 = r'*r;
			 beta = aux1/aux0;
			 p = r + beta*p;        % new search direction
			 aux0 = aux1;
         end
         
		% Step 2: Update the hyperparameter theta
		q = APChol'*spdiags(kron(sqrt(theta),ones(1,3))',0,3*N,3*N)*w;
		dip_norm2 = sum(reshape(q,3,N).^2,1);
		theta_old = theta;
		theta = 0.5*theta_star.*(eta + sqrt(eta^2 + 2*dip_norm2./theta_star));  % updated hyperparameter
		theta_diff=norm(theta_old-theta)/norm(theta_old);
		% Checking the theta condition
		if theta_diff<theta_tol 
		   % The relative change of theta is below the threshold; stop the outer iteration
		   outer_iteration = 'no';
        end
        diagnostics(outer_count,jj) = j;
        diagnostics(n_outer+1,jj) = theta_diff;
	 end
	 % Convergence is reached; saving dipole moments at time jj
	 Q(:,jj) = q;
end

