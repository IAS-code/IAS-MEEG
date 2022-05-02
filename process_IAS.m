function varargout = process_IAS( varargin )
% PROCESS_IAS: Call IAS algorithm
% 
% This function solves the MEEG inverse problem by the Iterative Alternating Sequential (IAS) algorithm. 
% The IAS algorithm is based on an iterative scheme that alternatively updates the dipole
% moments Q by solving a linear least squares problem using a priorconditioned CGLS algorithm with sutable 
% stopping condition and updating the hyperparameter theta by an explicit formula.
% 
% REFERENCES: https://ias-code.github.io/IAS-MEEG/index.html
% 
% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as publishebead by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: 

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'IAS';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 1000;
    sProcess.Description = 'https://ias-code.github.io/IAS-MEEG/index.html';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
    % Options: Time window
    sProcess.options.label1.Comment = '<B>Input options</B>:';
    sProcess.options.label1.Type    = 'label';
    sProcess.options.active_time.Comment = 'Time window: ';
    sProcess.options.active_time.Type    = 'timewindow';
    sProcess.options.active_time.Value   = [];
    % Options: Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, EEG';
    
    % Options: IAS parameters
    sProcess.options.label2.Comment = '<BR><B>IAS options</B>:';
    sProcess.options.label2.Type    = 'label';
    sProcess.options.SNR.Comment = 'bs_SNR (amplitude):';
    sProcess.options.SNR.Type    = 'value';
    sProcess.options.SNR.Value   = {3, ''};
    % Options: eta
    sProcess.options.eta.Comment = 'Focality parameter: ';
    sProcess.options.eta.Type    = 'text';
    sProcess.options.eta.Value   = '0.001';
    % Options: cut-off
    sProcess.options.cutoff.Comment = 'Cut off: ';
    sProcess.options.cutoff.Type    = 'value';
    sProcess.options.cutoff.Value   = {0.9, ''};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Initialize returned list of files
    OutputFiles = {};
    
    % ===== GET OPTIONS =====
    ActiveTime  = sProcess.options.active_time.Value{1}
    SensorTypes = sProcess.options.sensortypes.Value
    SNR         = sProcess.options.SNR.Value{1}^2
    Eta         = str2double(sProcess.options.eta.Value)
    CutOff      = sProcess.options.cutoff.Value{1}
%     
    
    % ===== LOAD CHANNEL FILE =====
    % Load channel file
    ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
    % Find the MEG channels
%     iMEG = good_channel(ChannelMat.Channel, [], 'MEG')
%     iEEG = good_channel(ChannelMat.Channel, [], 'EEG');
    iChannels = channel_find(ChannelMat.Channel, SensorTypes);
    
    % ===== LOAD HEAD MODEL =====
    % Get channel study
    [sChannelStudy, ~] = bst_get('ChannelFile', sInputs(1).ChannelFile);
    % Load the default head model
    HeadModelFile = sChannelStudy.HeadModel(sChannelStudy.iHeadModel).FileName;
    sHeadModel = load(file_fullpath(HeadModelFile));
    % Get number of sources
    nSources = length(sHeadModel.GridLoc);
    
    % ===== LOAD THE DATA =====
    DataMat = in_bst(sInputs(1).FileName, [], 0);
    nChannels = size(DataMat.F,1);
    nTime     = size(DataMat.F,2);
    Time      = DataMat.Time;
    
    iActiveTime = panel_time('GetTimeIndices', Time, ActiveTime);

    % Remove bad channels
    iBadChan = find(DataMat.ChannelFlag == -1);
    iChannels = setdiff(iChannels, iBadChan);
    % Error: All channels tagged as bad
    if isempty(iChannels)
        bst_report('Error', sProcess, sInput, 'All the selected channels are tagged as bad.');
        return;
    end
    
    % ===== COMPUTE ANATOMICAL PRIOR =====
    disp('Anatomical prior ...');
    APChol = BuildAnatomicalPrior(sHeadModel.GridLoc, sHeadModel.GridOrient);
    size(APChol)
    
    % ===== MODEL SCALING =====
    % Model scaling, definition of the scaling parameter theta^*
    disp('Model scaling ...');
    LF = sHeadModel.Gain(iChannels, :);
    B = DataMat.F(iChannels, iActiveTime);
    length(iChannels)
    size(LF)
    size(B)
    [theta_star, theta_cut_off, sigma, LF_scaling, B_scaling] = SetParameters(LF, APChol, B, SNR, CutOff);

    
    % ===== PROCESS =====
    % Run the IAS algorithm
    disp('Running IAS ...');
    [Q, diagnostics] = IAS_algorithm(LF, LF_scaling, APChol, B, B_scaling, sigma, theta_star, Eta);
    
    ImageGridAmp = zeros(3*nSources, nTime);
    ImageGridAmp(:, iActiveTime) = Q;

    % ===== SAVE THE RESULTS =====
    % Create a new data file structure
    ResultsMat = db_template('resultsmat');
    ResultsMat.ImagingKernel = [];
    ResultsMat.ImageGridAmp  = ImageGridAmp;
    ResultsMat.nComponents   = 3;   % 1 (constrained) or 3 (uncostrained)
    ResultsMat.Comment       = 'IAS';
    ResultsMat.Function      = 'IAS';
    ResultsMat.Time          = Time;  % Leave it empty if using ImagingKernel
    ResultsMat.iActiveTime   = 1:nTime; 
    ResultsMat.DataFile      = [];
    ResultsMat.HeadModelFile = HeadModelFile;
    ResultsMat.HeadModelType = sHeadModel.HeadModelType;
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = iChannels;
    ResultsMat.SurfaceFile   = sHeadModel.SurfaceFile;
    ResultsMat.GridLoc       = [];
    ResultsMat.SNR           = SNR;
    ResultsMat.CutOff        = CutOff;
    ResultsMat.Eta           = Eta;

    % === NOT SHARED ===
    % Get the output study (pick the one from the first file)
    iStudy = sInputs(1).iStudy;
    % Create a default output filename 
    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sInputs(1).FileName), 'results_');
%     % === SHARED ===
%     % Get the output study (pick the one from the first file)
%     iStudy = iChannelStudy;
%     % Create a default output filename 
%     OutputFiles{1} = bst_process('GetNewFilename', fileparts(sInputs(1).ChannelFile), 'results_KERNEL_');

%     OutputFiles{2} = ...

    % Save on disk
    save(OutputFiles{1}, '-struct', 'ResultsMat');
    % Register in database
    db_add_data(iStudy, OutputFiles{1}, ResultsMat);
end


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

  if size(coord, 1) ~= 3
    coord = coord';
  end
  if size(normals, 1) ~= 3
    normals = normals';
  end
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

end


function [theta_star,theta_cut_off,sigma,LF_scaling,B_scaling] = SetParameters(LF, APChol, B, SNR, cut_off, is_paint)
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
  [m, n3] = size(LF);
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
 
end


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

end

