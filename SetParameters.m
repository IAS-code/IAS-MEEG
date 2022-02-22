function [theta_star,theta_cut_off,sigma,LF_scaling,B_scaling] = SetParameters(LF,APChol,B,SNR,cut_off,is_paint)
% This function scales the lead field and the data, adjusts the truncation of the sensitivities, and returns an estimate for the standard deviation of the noise
%
%  Input:
%       LF: (m,3*n) array, the lead field matrix where m is the number of channels and n is the number of dipoles
%
%       APChol: Cholesky factor of the anatomical prior covariance       
%
%       B: (m,T) array, a clip of the data of length T
%
%       SNR: estimated signal-to-noise ratio
%
%       cut_off: scalar 0<cut_off<1, defines the cut off level, removing (1-cut_off)% of the highest theta_star values. Values of the order 0.8...0.9 are expected for cut_off (default value: 0.9)
%
%       is_paint: scalar; setting is_paint=1 plots theta_star (default value: 0)
%
%  Output:    
%       theta_star: n-vector, computed theta_star values with truncation
%
%       theta_cut_off: scalar, the cut-off value of the theta_star
%
%       sigma: scalar, standard deviation of the scaled noise
%
%       LF_scaling: scalar to be used to scale the lead field
%
%       B_scaling: scalar, to be used for scaling the data
%
%  Usage: 
%       [theta_star,theta_cut_off,sigma,LF_scaling,B_scaling] = SetParameters(LF,B,SNR,cut_off,is_paint)
%        SetParameters(LF,B,SNR) is equivalent to SetParameters(LF,B,SNR,0.9,0)
%
%        SetParameters(LF,B,SNR,cut_off) is equivalent to SetParameters(LF,B,SNR,cut_off,0)

% Setting default values
if nargin == 4, cut_off = 1; is_paint = 0; end
if nargin == 5,                is_paint = 0; end

% Reading the number of channels m and the number of dipoles n
[m,n3] = size(LF);
n = n3/3; 
%% 1. Lead field scaling

% Computing the sensitivity matrix S without anatomical prior: the jth row 
% contains the powers of unit dipoles at each field point, as registered by
% the jth magnetometer/electrode.
AddCols    = kron(speye(n),[1;1;1]);
S          = (LF.^2)*AddCols;
% Finding the power of the magnetic field/voltage at the most sensitive 
% magnetometer/electrode for each dipole
MaxS       = max(S);
% Scaling the lead field so that the mean maximum sensitivity is a unit
MeanSens   = 1/n*sum(MaxS);
% The lead field matrix is scaled by the square root (amplitude) of this number 
LF_scaling = 1/sqrt(MeanSens);
LF         = LF_scaling*LF;

%% 2. Scaling of the data

% Checking the size of the data
[m,T]     = size(B);
% Finding the average amplitude over the sensors for each time slice
Bnorm_ave = 1/m*sqrt(sum(B.^2,1))';
% Averaging over the time slices
MeanB     = 1/T*sum(Bnorm_ave);
% Scaling the data with the mean amplitude
B_scaling = 1/MeanB;
B         = B_scaling*B;

%% 3. Computing theta_star using the scaled fields

Phi        = 1/T*norm(B,'fro')^2;         % Mean power over time slices
S_scaled   = ((LF*APChol').^2)*AddCols;   % Scaled sensitivities
Sens       = sum(S_scaled,1);
% Safeguard: If SNR is given too low, theta_star ill-defined 
SNR_eff    = max(SNR,1.1);
theta_star = (2/5)*Phi*(1 - 1/SNR_eff)./Sens; % beta = 5/2

%% 4. Cutting possibly the highest theta_star values

Q             = theta_star - min(theta_star);
Qsort         = sort(Q,'ascend');
i_cut         = round(cut_off*n);
theta_cut_off = Qsort(i_cut) + min(theta_star);
theta_star    = min(theta_star,theta_cut_off);

% Optional: Plotting the values of theta_star
if is_paint
  figure
  semilogy((1:n),theta_star,'k.','MarkerSize',15)
  set(gca,'FontSize',20)
  ylabel('\theta^*,  normalized units','FontSize',20)
  % Adjusting the cut-off value
  hold on
  semilogy([1,n],[theta_cut_off,theta_cut_off],'r-','LineWidth',2)
  v = axis;
  axis([1,n,v(3),v(4)]);
  hold off
end

%% 5.  Adjusting the noise level of the scaled model

 B_ave_power = 1/T*sum(sum(B.^2)); 
 sigma2      = B_ave_power/(m*SNR);
 sigma       = sqrt(sigma2);
 
 