%% IAS_Demo
% This is a demo for testing the Iterative Alternating Sequential (IAS) algorithm. 
% 
% The real data set used for this demo is the MEG sample data acquired with the Neuromag Vectorview system 
% at MGH/HMS/MIT Athinoula A. Martinos Center Biomedical Imaging and made available, together with the 
% MRI reconstructions created with FreeSurfer, in the MNE software package (http://martinos.org/mne). 
% As part of the protocol for the data collection, checkerboard patterns were presented into the left 
% and right visual field, interspersed by tones to the left or right ear. The interval between the stimuli 
% was 750 ms. Occasionally a smiley face was presented at the center of the visual field. The subject 
% was asked to press a key with the right index finger as soon as possible after the appearance of the face 
% [A. Gramfort et al, MNE software for processing MEG and EEG data, Neuroimage 86 (2014) 446-460].
% In the test we only consider the trials corresponding to the left visual stimulus 
% and perform the averaging on these trials. 
% The source space containing both cortical regions and subcortical structures
% after discretization comprises 19054 vertices.
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

%% Reset all
clear; clc; 
% close all

%% Input
%     SNR: estimated signal-to-noise ratio
SNR = 9;
%     eta: scalar, parameter for selecting the focality of the reconstructed activity
eta = 0.01;
cut_off = 0.9;
%% Loading the source space
%     coord, normals: (3,N) array, coordinates and normal vectors of the dipoles in the grid
disp('Loading source space')
load('SourceSpace_DBA')

%% Loading the leadfield matrix
%     LF: (M,3*N) array, the lead field matrix (M is the number of channels; N is the number of dipoles)
disp('Loading leadfield matrix')
load('LeadfieldMatrix_DBA')

%% Loading the magnetic data
%     data: (M,T) array, a set of data of length T
%     time: T-vector, times 
disp('Loading data')
load('visualMEGData')

%% Building anatomical prior
disp('Building anatomical prior')
APChol = BuildAnatomicalPrior(coord,normals);

%% Setting theta_star and scaling 
disp('Setting parameters')
t_peak = 86;   
t_min = 1;
t_max = 211;
B = data(:,t_min:t_max);
[theta_star,theta_cut_off,sigma,LF_scaling,B_scaling] = SetParameters(LF,APChol,B,SNR,cut_off);

%% Solving the inverse problem using IAS algorithm
disp(['Running IAS algorithm from time ',num2str(time(t_min)),' ms to time ',num2str(time(t_max)),' ms'])
Q = IAS_algorithm(LF, LF_scaling, APChol, B, B_scaling, sigma, theta_star, eta);

%% Visualizing the activity map
disp(['Visualizing activity map at time ',num2str(time(t_peak)),' ms'])
N = size(LF,2)/3;  % number of dipoles in the source space
t_vis = t_peak - t_min +1;
q = Q(:, t_vis);
dip_norm2 = sum(reshape(q,3,N).^2,1);
Q_est = sqrt(dip_norm2);
% SlicedVisualization_ActivityMap(coord,Q_est);


%% Visualazing the activity map and its maximum
[Q_max,i_max] = max(Q_est);
r0 = coord(:,i_max);
SingleSliceVisualization_ActivityMap(coord,Q_est,r0,'CoordSystem','Brainstorm','ColorScale','double','Locator','pinhead');


%% Compute time series of maximum
% Time series of all dipoles of the source space
nt = t_max - t_min + 1;

Q_est_all = zeros(N, nt);
for i=1:size(Q,2)
  q = Q(:, i);
  Q_est_all(:, i) = sqrt(sum(reshape(q,3,N).^2,1));
end
figure, plot(Q_est_all')


%% Visualizing time series
t_sel = [83, 110, 128, 143, 201] - t_min +1;  % visual selection of tp
color = {[255 140   0]/255, 'r', [139   0   0]/255, 'c', [0   0 255]/255, [30 144 255]/255, [128   0 128]/255};

figure
subplot(1,2,1)
plot3(coord(1, :), coord(2, :), coord(3, :),'color', [128 128 128]/255, 'Marker', '.', 'MarkerSize', 4, 'linestyle', 'none')
hold on
subplot(1, 2, 2), hold on
for i=1:length(t_sel)
  t = t_sel(i);
  time(t_sel(i))
  [val1, i_val1] = max(Q_est_all(:, t));
  
  subplot(1,2,1), plot3(coord(1, i_val1), coord(2, i_val1), coord(3, i_val1), 'color', color{i} , 'Marker', '.', 'MarkerSize', 45)
  subplot(1, 2, 2), plot(time(t_min:t_max), Q_est_all(i_val1, :), 'color', color{i}, 'LineWidth', 2)
  
end

print('time_series.png', '-dpng', '-r300')


