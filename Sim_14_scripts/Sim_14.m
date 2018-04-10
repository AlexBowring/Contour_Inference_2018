function Sim_14(nSubj,SvNm,nRlz)
%
% Creates a 2D images of linearly increasing signal from L to R, and then applies the standardized effects Contour Inference method
% for each of the proposed options
%


%------------Starting Up initialization
if (nargin<1)
  nSubj  = 60;  % Number of subjects
end
if (nargin<2)
  SvNm  = 'Normsim';  % Save name
end
if (nargin<3)
  nRlz = 5000;
end  
if exist([SvNm '.mat'], 'file')
  error('Will not overwrite sim result')
end

%------------Define parameters
% SvNm = 'LinearSig';
% nSubj  = 120;
% nRlz = 300;

tau     = 1/sqrt(nSubj);
nBoot   = 5000;
dim     = [100 100]; 
mag     = 3;
smo     = 10;
rimFWHM = 15; 				 
stdblk  = prod(dim([1 2])/2);
thr     = 2;

% Adding variable for low resolutions simulation
low_res_dim = [50 50];

%-----------Initialization of Some Variables
V           = prod(dim);   
wdim        = dim + 2*ceil(rimFWHM*smo*ones(1,2));  % Working image dimension
trunc_x     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(1))};
trunc_y     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(2))};
trnind      = cat(2, trunc_x, trunc_y);

observed_data  = zeros([dim nSubj]);
raw_noise      = zeros([wdim nSubj]);

% This stores the vector SupG for each run
% This vector stores the result for each realisation on whether AC^+ < AC < AC^ for each level of smoothing (1 if true, 0 if false) 
subset_success_vector_raw_80           = zeros(nRlz, 1); 
subset_success_vector_raw_90           = zeros(nRlz, 1);
subset_success_vector_raw_95           = zeros(nRlz, 1);

subset_success_vector_raw_80_linear    = zeros(nRlz, 1); 
subset_success_vector_raw_90_linear    = zeros(nRlz, 1);
subset_success_vector_raw_95_linear    = zeros(nRlz, 1);

low_res_subset_success_vector_raw_80           = zeros(nRlz, 1); 
low_res_subset_success_vector_raw_90           = zeros(nRlz, 1);
low_res_subset_success_vector_raw_95           = zeros(nRlz, 1);

low_res_subset_success_vector_raw_80_linear    = zeros(nRlz, 1); 
low_res_subset_success_vector_raw_90_linear    = zeros(nRlz, 1);
low_res_subset_success_vector_raw_95_linear    = zeros(nRlz, 1); 

%- This vector stores the threshold value 'c' for each run
threshold_raw_80_store                  = zeros(nRlz, 1);
threshold_raw_90_store                  = zeros(nRlz, 1);
threshold_raw_95_store                  = zeros(nRlz, 1);

threshold_raw_80_linear_store           = zeros(nRlz, 1);
threshold_raw_90_linear_store           = zeros(nRlz, 1);
threshold_raw_95_linear_store           = zeros(nRlz, 1);

low_res_threshold_raw_80_store                  = zeros(nRlz, 1);
low_res_threshold_raw_90_store                  = zeros(nRlz, 1);
low_res_threshold_raw_95_store                  = zeros(nRlz, 1);

low_res_threshold_raw_80_linear_store           = zeros(nRlz, 1);
low_res_threshold_raw_90_linear_store           = zeros(nRlz, 1);
low_res_threshold_raw_95_linear_store           = zeros(nRlz, 1);

%- This vector stores the percentage volumes A^+_c/A_c, A^_c/A_c, A^-_c/A_c
lower_contour_raw_80_volume_prct_store          = zeros(nRlz, 1);
upper_contour_raw_80_volume_prct_store          = zeros(nRlz, 1);
lower_contour_raw_80_linear_volume_prct_store   = zeros(nRlz, 1);
upper_contour_raw_80_linear_volume_prct_store   = zeros(nRlz, 1);

lower_contour_raw_90_volume_prct_store          = zeros(nRlz, 1);
upper_contour_raw_90_volume_prct_store          = zeros(nRlz, 1);
lower_contour_raw_90_linear_volume_prct_store   = zeros(nRlz, 1);
upper_contour_raw_90_linear_volume_prct_store   = zeros(nRlz, 1);

lower_contour_raw_95_volume_prct_store          = zeros(nRlz, 1);
upper_contour_raw_95_volume_prct_store          = zeros(nRlz, 1);
lower_contour_raw_95_linear_volume_prct_store   = zeros(nRlz, 1);
upper_contour_raw_95_linear_volume_prct_store   = zeros(nRlz, 1);

low_res_lower_contour_raw_80_volume_prct_store          = zeros(nRlz, 1);
low_res_upper_contour_raw_80_volume_prct_store          = zeros(nRlz, 1);
low_res_lower_contour_raw_80_linear_volume_prct_store   = zeros(nRlz, 1);
low_res_upper_contour_raw_80_linear_volume_prct_store   = zeros(nRlz, 1);

low_res_lower_contour_raw_90_volume_prct_store          = zeros(nRlz, 1);
low_res_upper_contour_raw_90_volume_prct_store          = zeros(nRlz, 1);
low_res_lower_contour_raw_90_linear_volume_prct_store   = zeros(nRlz, 1);
low_res_upper_contour_raw_90_linear_volume_prct_store   = zeros(nRlz, 1);

low_res_lower_contour_raw_95_volume_prct_store          = zeros(nRlz, 1);
low_res_upper_contour_raw_95_volume_prct_store          = zeros(nRlz, 1);
low_res_lower_contour_raw_95_linear_volume_prct_store   = zeros(nRlz, 1);
low_res_upper_contour_raw_95_linear_volume_prct_store   = zeros(nRlz, 1);

% This stores the vector SupG for each run
supG_raw_store         = zeros(nBoot, nRlz);
supG_raw_linear_store  = zeros(nBoot, nRlz);
low_res_supG_raw_store         = zeros(nBoot, nRlz);
low_res_supG_raw_linear_store  = zeros(nBoot, nRlz);

%-These matrices store all the sets of interest during the bootstrap
% method for all levels of smoothing
lower_contour_raw_80_store                       = zeros([nRlz dim]);
upper_contour_raw_80_store                       = zeros([nRlz dim]);
upper_subset_mid_raw_80_store                    = zeros([nRlz dim]);
mid_subset_lower_raw_80_store                    = zeros([nRlz dim]);
lower_contour_raw_80_linear_store                = zeros([nRlz dim]);
upper_contour_raw_80_linear_store                = zeros([nRlz dim]);
upper_subset_mid_raw_80_linear_store             = zeros([nRlz dim]);
mid_subset_lower_raw_80_linear_store             = zeros([nRlz dim]);

low_res_lower_contour_raw_80_store                       = zeros([nRlz low_res_dim]);
low_res_upper_contour_raw_80_store                       = zeros([nRlz low_res_dim]);
low_res_upper_subset_mid_raw_80_store                    = zeros([nRlz low_res_dim]);
low_res_mid_subset_lower_raw_80_store                    = zeros([nRlz low_res_dim]);
low_res_lower_contour_raw_80_linear_store                = zeros([nRlz low_res_dim]);
low_res_upper_contour_raw_80_linear_store                = zeros([nRlz low_res_dim]);
low_res_upper_subset_mid_raw_80_linear_store             = zeros([nRlz low_res_dim]);
low_res_mid_subset_lower_raw_80_linear_store             = zeros([nRlz low_res_dim]);

lower_contour_raw_90_store                       = zeros([nRlz dim]);
upper_contour_raw_90_store                       = zeros([nRlz dim]);
upper_subset_mid_raw_90_store                    = zeros([nRlz dim]);
mid_subset_lower_raw_90_store                    = zeros([nRlz dim]);
lower_contour_raw_90_linear_store                = zeros([nRlz dim]);
upper_contour_raw_90_linear_store                = zeros([nRlz dim]);
upper_subset_mid_raw_90_linear_store             = zeros([nRlz dim]);
mid_subset_lower_raw_90_linear_store             = zeros([nRlz dim]);

low_res_lower_contour_raw_90_store                       = zeros([nRlz low_res_dim]);
low_res_upper_contour_raw_90_store                       = zeros([nRlz low_res_dim]);
low_res_upper_subset_mid_raw_90_store                    = zeros([nRlz low_res_dim]);
low_res_mid_subset_lower_raw_90_store                    = zeros([nRlz low_res_dim]);
low_res_lower_contour_raw_90_linear_store                = zeros([nRlz low_res_dim]);
low_res_upper_contour_raw_90_linear_store                = zeros([nRlz low_res_dim]);
low_res_upper_subset_mid_raw_90_linear_store             = zeros([nRlz low_res_dim]);
low_res_mid_subset_lower_raw_90_linear_store             = zeros([nRlz low_res_dim]);

lower_contour_raw_95_store                       = zeros([nRlz dim]);
upper_contour_raw_95_store                       = zeros([nRlz dim]);
upper_subset_mid_raw_95_store                    = zeros([nRlz dim]);
mid_subset_lower_raw_95_store                    = zeros([nRlz dim]);
lower_contour_raw_95_linear_store                = zeros([nRlz dim]);
upper_contour_raw_95_linear_store                = zeros([nRlz dim]);
upper_subset_mid_raw_95_linear_store             = zeros([nRlz dim]);
mid_subset_lower_raw_95_linear_store             = zeros([nRlz dim]);

low_res_lower_contour_raw_95_store                       = zeros([nRlz low_res_dim]);
low_res_upper_contour_raw_95_store                       = zeros([nRlz low_res_dim]);
low_res_upper_subset_mid_raw_95_store                    = zeros([nRlz low_res_dim]);
low_res_mid_subset_lower_raw_95_store                    = zeros([nRlz low_res_dim]);
low_res_lower_contour_raw_95_linear_store                = zeros([nRlz low_res_dim]);
low_res_upper_contour_raw_95_linear_store                = zeros([nRlz low_res_dim]);
low_res_upper_subset_mid_raw_95_linear_store             = zeros([nRlz low_res_dim]);
low_res_mid_subset_lower_raw_95_linear_store             = zeros([nRlz low_res_dim]);

supG_raw              = zeros(nBoot,1);
supG_raw_linear       = zeros(nBoot,1);
low_res_supG_raw              = zeros(nBoot,1);
low_res_supG_raw_linear       = zeros(nBoot,1);

% Creating linearly increasing signal across columns
Sig = repmat(linspace(1, 3), dim(2), 1);

% Uncomment to look at the Signal
%imagesc(Sig); axis image; colorbar
AC = Sig >= thr;

% Variables for computing the estimated boundary
[a,b] = ndgrid(-1:1);
se = strel('arbitrary',sqrt(a.^2 + b.^2) <=1);

% The boundary for Sig > 2, note that Sig = 2.02 in the 51st column
true_boundary = zeros(dim);
true_boundary(:,51) = ones(100, 1);
true_boundary = logical(true_boundary);

% low resolution signal and boundary
low_res_Sig                 = Sig(2:2:end,2:2:end);
low_res_AC                  = low_res_Sig >= thr; 
low_res_true_boundary       = zeros(low_res_dim);
low_res_true_boundary(:,26) = ones(50, 1);
low_res_true_boundary       = logical(low_res_true_boundary);

for t=1:nRlz
    fprintf('.');
      for i=1:nSubj
	    %
	    % Generate random realizations of signal + noise
	    %
        raw_noise(:,:,i) = randn(wdim); %- Noise that will be added to the signal 

        %
        % smooth noise  
        %
        [Noises,tt] = spm_conv(raw_noise(:,:,i),smo,smo);
        Noises = Noises/sqrt(tt);      
      
        %
        % Truncate to avoid edge effects
        %
        tNoises = Noises(trnind{1},trnind{2});       
        tImgs = Sig + tNoises; % Creates the true image of smoothed signal + smoothed noise
        observed_data(:,:,i) = tImgs;
        
      end %========== Loop i (subjects)
      
      
      %% High res variables
      observed_mean = mean(observed_data,3);

      observed_std = reshape(...
         biasmystd(reshape(observed_data,[prod(dim) nSubj]),stdblk),...
           dim);
       
      % Making the three observed boundaries: dilated boundary, eroded
      % boundary, and dilated - eroded boundary.
      observed_AC = observed_mean >= thr;
      observed_AC_volume = sum(observed_AC(:)); 

      % Making the interpolated boundary edges
      % Horizontal edges
      horz = observed_AC(:,2:end) | observed_AC(:,1:end-1);
      % Compute the left shifted horizontal edges
      lshift            = observed_AC; % initialize
      lshift(:,1:end-1) = horz;
      lshift            = lshift & ~observed_AC;
      %%% Compute the right shifted horizontal edges
      rshift          = observed_AC; % initialize
      rshift(:,2:end) = horz;
      rshift          = rshift & ~observed_AC;
      % Vertical edges
      vert = observed_AC(1:end-1,:) | observed_AC(2:end,:);
      %%% Compute the right shifted horizontal edges
      ushift = observed_AC;
      ushift(1:end-1,:) = vert;
      ushift = ushift & ~observed_AC;
      %%% Compute the down shifted vertical edges
      dshift = observed_AC;
      %%% Values of random field on down shifted vertical edges
      dshift(2:end,:)   = vert;
      dshift = dshift & ~observed_AC;
 
      % Residuals
      resid = bsxfun(@minus,observed_data,observed_mean);
      resid = spdiags(1./reshape(observed_std, [prod(dim) 1]), 0,prod(dim),prod(dim))*reshape(resid,[prod(dim) nSubj]); 
       
      %% Low res variables
      low_res_observed_data = observed_data(2:2:end,2:2:end,:);
      
            low_res_observed_mean = mean(low_res_observed_data,3);

      low_res_observed_std = reshape(...
         biasmystd(reshape(low_res_observed_data,[prod(low_res_dim) nSubj]),stdblk),...
           low_res_dim);
       
      % Making the three observed boundaries: dilated boundary, eroded
      % boundary, and dilated - eroded boundary.
      low_res_observed_AC = low_res_observed_mean >= thr;
      low_res_observed_AC_volume = sum(low_res_observed_AC(:)); 

      % Making the interpolated boundary edges
      % Horizontal edges
      low_res_horz = low_res_observed_AC(:,2:end) | low_res_observed_AC(:,1:end-1);
      % Compute the left shifted horizontal edges
      low_res_lshift            = low_res_observed_AC; % initialize
      low_res_lshift(:,1:end-1) = low_res_horz;
      low_res_lshift            = low_res_lshift & ~low_res_observed_AC;
      %%% Compute the right shifted horizontal edges
      low_res_rshift          = low_res_observed_AC; % initialize
      low_res_rshift(:,2:end) = low_res_horz;
      low_res_rshift          = low_res_rshift & ~low_res_observed_AC;
      % Vertical edges
      low_res_vert = low_res_observed_AC(1:end-1,:) | low_res_observed_AC(2:end,:);
      %%% Compute the right shifted horizontal edges
      low_res_ushift = low_res_observed_AC;
      low_res_ushift(1:end-1,:) = low_res_vert;
      low_res_ushift = low_res_ushift & ~low_res_observed_AC;
      %%% Compute the down shifted vertical edges
      low_res_dshift = low_res_observed_AC;
      %%% Values of random field on down shifted vertical edges
      low_res_dshift(2:end,:)   = low_res_vert;
      low_res_dshift = low_res_dshift & ~low_res_observed_AC;
 
      % Residuals
      low_res_resid = bsxfun(@minus,low_res_observed_data,low_res_observed_mean);
      low_res_resid = spdiags(1./reshape(low_res_observed_std, [prod(low_res_dim) 1]), 0,prod(low_res_dim),prod(low_res_dim))*reshape(low_res_resid,[prod(low_res_dim) nSubj]); 
      
      %% Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          % Applying the bootstrap using Rademacher variables (signflips)
          signflips                              = randi(2,[nSubj,1])*2-3;
          
          %% High res case
          resid_bootstrap                        = resid*spdiags(signflips, 0, nSubj, nSubj);
          resid_bootstrap                        = reshape(resid_bootstrap, [dim nSubj]);
          resid_field                            = sum(resid_bootstrap, 3)/sqrt(nSubj); 

          supG_raw(k)          = max(abs(resid_field(true_boundary)));
          
          % Calculating the maximum over the linear boundary edges
          lshift_boundary_values = abs((resid_field(lshift) + resid_field(lshift(:,[dim(2) 1:dim(2)-1])))/2);
          rshift_boundary_values = abs((resid_field(rshift) + resid_field(rshift(:,[2:dim(2) 1])))/2);
          ushift_boundary_values = abs((resid_field(ushift) + resid_field(ushift([dim(1) 1:dim(1)-1],:)))/2);
          dshift_boundary_values = abs((resid_field(dshift) + resid_field(dshift([2:dim(1) 1],:)))/2);
          supG_raw_linear(k)   = max([lshift_boundary_values; rshift_boundary_values; ushift_boundary_values; dshift_boundary_values]);
          
          %% Low res case
          low_res_resid_bootstrap                        = low_res_resid*spdiags(signflips, 0, nSubj, nSubj);
          low_res_resid_bootstrap                        = reshape(low_res_resid_bootstrap, [low_res_dim nSubj]);
          low_res_resid_field                            = sum(low_res_resid_bootstrap, 3)/sqrt(nSubj); 

          low_res_supG_raw(k)          = max(abs(low_res_resid_field(low_res_true_boundary)));
          
          % Calculating the maximum over the linear boundary edges
          low_res_lshift_boundary_values = abs((low_res_resid_field(low_res_lshift) + low_res_resid_field(low_res_lshift(:,[low_res_dim(2) 1:low_res_dim(2)-1])))/2);
          low_res_rshift_boundary_values = abs((low_res_resid_field(low_res_rshift) + low_res_resid_field(low_res_rshift(:,[2:low_res_dim(2) 1])))/2);
          low_res_ushift_boundary_values = abs((low_res_resid_field(low_res_ushift) + low_res_resid_field(low_res_ushift([low_res_dim(1) 1:low_res_dim(1)-1],:)))/2);
          low_res_dshift_boundary_values = abs((low_res_resid_field(low_res_dshift) + low_res_resid_field(low_res_dshift([2:low_res_dim(1) 1],:)))/2);
          low_res_supG_raw_linear(k)   = max([low_res_lshift_boundary_values; low_res_rshift_boundary_values; low_res_ushift_boundary_values; low_res_dshift_boundary_values]);
      end
    
    %% High res case
    middle_contour                = AC;
    middle_contour_volume         = sum(middle_contour(:));
    
    % Gaussian random variable results for the true and estimated boundary
    % True boundary
    supGa_raw_80                     = prctile(supG_raw, 80);
    supGa_raw_90                     = prctile(supG_raw, 90);
    supGa_raw_95                     = prctile(supG_raw, 95);
       
    lower_contour_raw_80             = observed_mean >= thr - supGa_raw_80*tau*observed_std;
    upper_contour_raw_80             = observed_mean >= thr + supGa_raw_80*tau*observed_std;
    lower_contour_raw_80_volume_prct = sum(lower_contour_raw_80(:))/middle_contour_volume;
    upper_contour_raw_80_volume_prct = sum(upper_contour_raw_80(:))/middle_contour_volume;
    mid_on_upper_raw_80              = upper_contour_raw_80.*middle_contour;
    lower_on_mid_raw_80              = middle_contour.*lower_contour_raw_80;
    upper_subset_mid_raw_80          = upper_contour_raw_80 - mid_on_upper_raw_80;
    mid_subset_lower_raw_80          = middle_contour - lower_on_mid_raw_80;
    
    lower_contour_raw_90             = observed_mean >= thr - supGa_raw_90*tau*observed_std;
    upper_contour_raw_90             = observed_mean >= thr + supGa_raw_90*tau*observed_std;
    lower_contour_raw_90_volume_prct = sum(lower_contour_raw_90(:))/middle_contour_volume;
    upper_contour_raw_90_volume_prct = sum(upper_contour_raw_90(:))/middle_contour_volume;
    mid_on_upper_raw_90              = upper_contour_raw_90.*middle_contour;
    lower_on_mid_raw_90              = middle_contour.*lower_contour_raw_90;
    upper_subset_mid_raw_90          = upper_contour_raw_90 - mid_on_upper_raw_90;
    mid_subset_lower_raw_90          = middle_contour - lower_on_mid_raw_90;    
    
    lower_contour_raw_95             = observed_mean >= thr - supGa_raw_95*tau*observed_std;
    upper_contour_raw_95             = observed_mean >= thr + supGa_raw_95*tau*observed_std;
    lower_contour_raw_95_volume_prct = sum(lower_contour_raw_95(:))/middle_contour_volume;
    upper_contour_raw_95_volume_prct = sum(upper_contour_raw_95(:))/middle_contour_volume;
    mid_on_upper_raw_95              = upper_contour_raw_95.*middle_contour;
    lower_on_mid_raw_95              = middle_contour.*lower_contour_raw_95;
    upper_subset_mid_raw_95          = upper_contour_raw_95 - mid_on_upper_raw_95;
    mid_subset_lower_raw_95          = middle_contour - lower_on_mid_raw_95;

    % Linear Boundary
    supGa_raw_80_linear                     = prctile(supG_raw_linear, 80);
    supGa_raw_90_linear                     = prctile(supG_raw_linear, 90);
    supGa_raw_95_linear                     = prctile(supG_raw_linear, 95);
       
    lower_contour_raw_80_linear             = observed_mean >= thr - supGa_raw_80_linear*tau*observed_std;
    upper_contour_raw_80_linear             = observed_mean >= thr + supGa_raw_80_linear*tau*observed_std;
    lower_contour_raw_80_linear_volume_prct = sum(lower_contour_raw_80_linear(:))/middle_contour_volume;
    upper_contour_raw_80_linear_volume_prct = sum(upper_contour_raw_80_linear(:))/middle_contour_volume;
    mid_on_upper_raw_80_linear              = upper_contour_raw_80_linear.*middle_contour;
    lower_on_mid_raw_80_linear              = middle_contour.*lower_contour_raw_80_linear;
    upper_subset_mid_raw_80_linear          = upper_contour_raw_80_linear - mid_on_upper_raw_80_linear;
    mid_subset_lower_raw_80_linear          = middle_contour - lower_on_mid_raw_80_linear;
    
    lower_contour_raw_90_linear             = observed_mean >= thr - supGa_raw_90_linear*tau*observed_std;
    upper_contour_raw_90_linear             = observed_mean >= thr + supGa_raw_90_linear*tau*observed_std;
    lower_contour_raw_90_linear_volume_prct = sum(lower_contour_raw_90_linear(:))/middle_contour_volume;
    upper_contour_raw_90_linear_volume_prct = sum(upper_contour_raw_90_linear(:))/middle_contour_volume;
    mid_on_upper_raw_90_linear              = upper_contour_raw_90_linear.*middle_contour;
    lower_on_mid_raw_90_linear              = middle_contour.*lower_contour_raw_90_linear;
    upper_subset_mid_raw_90_linear          = upper_contour_raw_90_linear - mid_on_upper_raw_90_linear;
    mid_subset_lower_raw_90_linear          = middle_contour - lower_on_mid_raw_90_linear;    
    
    lower_contour_raw_95_linear             = observed_mean >= thr - supGa_raw_95_linear*tau*observed_std;
    upper_contour_raw_95_linear             = observed_mean >= thr + supGa_raw_95_linear*tau*observed_std;
    lower_contour_raw_95_linear_volume_prct = sum(lower_contour_raw_95_linear(:))/middle_contour_volume;
    upper_contour_raw_95_linear_volume_prct = sum(upper_contour_raw_95_linear(:))/middle_contour_volume;
    mid_on_upper_raw_95_linear              = upper_contour_raw_95_linear.*middle_contour;
    lower_on_mid_raw_95_linear              = middle_contour.*lower_contour_raw_95_linear;
    upper_subset_mid_raw_95_linear          = upper_contour_raw_95_linear - mid_on_upper_raw_95_linear;
    mid_subset_lower_raw_95_linear          = middle_contour - lower_on_mid_raw_95_linear;

    %
    % Storing all variables of interest
    %
    supG_raw_store(:,t)                                    = supG_raw;
    threshold_raw_80_store(t)                              = supGa_raw_80;
    lower_contour_raw_80_store(t,:,:)                      = lower_contour_raw_80;
    upper_contour_raw_80_store(t,:,:)                      = upper_contour_raw_80;
    upper_subset_mid_raw_80_store(t,:,:)                   = upper_subset_mid_raw_80;
    mid_subset_lower_raw_80_store(t,:,:)                   = mid_subset_lower_raw_80;
    lower_contour_raw_80_volume_prct_store(t)              = lower_contour_raw_80_volume_prct;
    upper_contour_raw_80_volume_prct_store(t)              = upper_contour_raw_80_volume_prct;
 
    threshold_raw_90_store(t)                              = supGa_raw_90;
    lower_contour_raw_90_store(t,:,:)                      = lower_contour_raw_90;
    upper_contour_raw_90_store(t,:,:)                      = upper_contour_raw_90;
    upper_subset_mid_raw_90_store(t,:,:)                   = upper_subset_mid_raw_90;
    mid_subset_lower_raw_90_store(t,:,:)                   = mid_subset_lower_raw_90;
    lower_contour_raw_90_volume_prct_store(t)              = lower_contour_raw_90_volume_prct;
    upper_contour_raw_90_volume_prct_store(t)              = upper_contour_raw_90_volume_prct;

    threshold_raw_95_store(t)                              = supGa_raw_95;
    lower_contour_raw_95_store(t,:,:)                      = lower_contour_raw_95;
    upper_contour_raw_95_store(t,:,:)                      = upper_contour_raw_95;
    upper_subset_mid_raw_95_store(t,:,:)                   = upper_subset_mid_raw_95;
    mid_subset_lower_raw_95_store(t,:,:)                   = mid_subset_lower_raw_95;
    lower_contour_raw_95_volume_prct_store(t)              = lower_contour_raw_95_volume_prct;
    upper_contour_raw_95_volume_prct_store(t)              = upper_contour_raw_95_volume_prct;

    supG_raw_linear_store(:,t)                                    = supG_raw_linear;
    threshold_raw_80_linear_store(t)                              = supGa_raw_80_linear;
    lower_contour_raw_80_linear_store(t,:,:)                      = lower_contour_raw_80_linear;
    upper_contour_raw_80_linear_store(t,:,:)                      = upper_contour_raw_80_linear;
    upper_subset_mid_raw_80_linear_store(t,:,:)                   = upper_subset_mid_raw_80_linear;
    mid_subset_lower_raw_80_linear_store(t,:,:)                   = mid_subset_lower_raw_80_linear;
    lower_contour_raw_80_linear_volume_prct_store(t)              = lower_contour_raw_80_linear_volume_prct;
    upper_contour_raw_80_linear_volume_prct_store(t)              = upper_contour_raw_80_linear_volume_prct;
 
    threshold_raw_90_linear_store(t)                              = supGa_raw_90_linear;
    lower_contour_raw_90_linear_store(t,:,:)                      = lower_contour_raw_90_linear;
    upper_contour_raw_90_linear_store(t,:,:)                      = upper_contour_raw_90_linear;
    upper_subset_mid_raw_90_linear_store(t,:,:)                   = upper_subset_mid_raw_90_linear;
    mid_subset_lower_raw_90_linear_store(t,:,:)                   = mid_subset_lower_raw_90_linear;
    lower_contour_raw_90_linear_volume_prct_store(t)              = lower_contour_raw_90_linear_volume_prct;
    upper_contour_raw_90_linear_volume_prct_store(t)              = upper_contour_raw_90_linear_volume_prct;

    threshold_raw_95_linear_store(t)                              = supGa_raw_95_linear;
    lower_contour_raw_95_linear_store(t,:,:)                      = lower_contour_raw_95_linear;
    upper_contour_raw_95_linear_store(t,:,:)                      = upper_contour_raw_95_linear;
    upper_subset_mid_raw_95_linear_store(t,:,:)                   = upper_subset_mid_raw_95_linear;
    mid_subset_lower_raw_95_linear_store(t,:,:)                   = mid_subset_lower_raw_95_linear;
    lower_contour_raw_95_linear_volume_prct_store(t)              = lower_contour_raw_95_linear_volume_prct;
    upper_contour_raw_95_linear_volume_prct_store(t)              = upper_contour_raw_95_linear_volume_prct;
    
    
    if sum(upper_subset_mid_raw_80(:))+sum(mid_subset_lower_raw_80(:))==0
      subset_success_vector_raw_80(t) = 1; 
      fprintf('raw nominal 80 true boundary success! \n');
    else 
      subset_success_vector_raw_80(t) = 0; 
      fprintf('raw nominal 80 true boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_90(:))+sum(mid_subset_lower_raw_90(:))==0
      subset_success_vector_raw_90(t) = 1; 
      fprintf('raw nominal 90 true boundary success! \n');
    else 
      subset_success_vector_raw_90(t) = 0; 
      fprintf('raw nominal 90 true boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_95(:))+sum(mid_subset_lower_raw_95(:))==0
      subset_success_vector_raw_95(t) = 1; 
      fprintf('raw nominal 95 true boundary success! \n');
    else 
      subset_success_vector_raw_95(t) = 0; 
      fprintf('raw nominal 95 true boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_80_linear(:))+sum(mid_subset_lower_raw_80_linear(:))==0
      subset_success_vector_raw_80_linear(t) = 1; 
      fprintf('raw nominal 80 linear boundary success! \n');
    else 
      subset_success_vector_raw_80_linear(t) = 0; 
      fprintf('raw nominal 80 linear boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_90_linear(:))+sum(mid_subset_lower_raw_90_linear(:))==0
      subset_success_vector_raw_90_linear(t) = 1; 
      fprintf('raw nominal 90 linear boundary success! \n');
    else 
      subset_success_vector_raw_90_linear(t) = 0; 
      fprintf('raw nominal 90 linear boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_95_linear(:))+sum(mid_subset_lower_raw_95_linear(:))==0
      subset_success_vector_raw_95_linear(t) = 1; 
      fprintf('raw nominal 95 linear boundary success! \n');
    else 
      subset_success_vector_raw_95_linear(t) = 0; 
      fprintf('raw nominal 95 linear boundary failure! \n');
    end     

    %% Low res case
    low_res_middle_contour                = low_res_AC;
    low_res_middle_contour_volume         = sum(low_res_middle_contour(:));
    
    % Gaussian random variable results for the true and estimated boundary
    % True boundary
    low_res_supGa_raw_80                     = prctile(low_res_supG_raw, 80);
    low_res_supGa_raw_90                     = prctile(low_res_supG_raw, 90);
    low_res_supGa_raw_95                     = prctile(low_res_supG_raw, 95);
       
    low_res_lower_contour_raw_80             = low_res_observed_mean >= thr - low_res_supGa_raw_80*tau*low_res_observed_std;
    low_res_upper_contour_raw_80             = low_res_observed_mean >= thr + low_res_supGa_raw_80*tau*low_res_observed_std;
    low_res_lower_contour_raw_80_volume_prct = sum(low_res_lower_contour_raw_80(:))/low_res_middle_contour_volume;
    low_res_upper_contour_raw_80_volume_prct = sum(low_res_upper_contour_raw_80(:))/low_res_middle_contour_volume;
    low_res_mid_on_upper_raw_80              = low_res_upper_contour_raw_80.*low_res_middle_contour;
    low_res_lower_on_mid_raw_80              = low_res_middle_contour.*low_res_lower_contour_raw_80;
    low_res_upper_subset_mid_raw_80          = low_res_upper_contour_raw_80 - low_res_mid_on_upper_raw_80;
    low_res_mid_subset_lower_raw_80          = low_res_middle_contour - low_res_lower_on_mid_raw_80;
    
    low_res_lower_contour_raw_90             = low_res_observed_mean >= thr - low_res_supGa_raw_90*tau*low_res_observed_std;
    low_res_upper_contour_raw_90             = low_res_observed_mean >= thr + low_res_supGa_raw_90*tau*low_res_observed_std;
    low_res_lower_contour_raw_90_volume_prct = sum(low_res_lower_contour_raw_90(:))/low_res_middle_contour_volume;
    low_res_upper_contour_raw_90_volume_prct = sum(low_res_upper_contour_raw_90(:))/low_res_middle_contour_volume;
    low_res_mid_on_upper_raw_90              = low_res_upper_contour_raw_90.*low_res_middle_contour;
    low_res_lower_on_mid_raw_90              = low_res_middle_contour.*low_res_lower_contour_raw_90;
    low_res_upper_subset_mid_raw_90          = low_res_upper_contour_raw_90 - low_res_mid_on_upper_raw_90;
    low_res_mid_subset_lower_raw_90          = low_res_middle_contour - low_res_lower_on_mid_raw_90;  
    
    low_res_lower_contour_raw_95             = low_res_observed_mean >= thr - low_res_supGa_raw_95*tau*low_res_observed_std;
    low_res_upper_contour_raw_95             = low_res_observed_mean >= thr + low_res_supGa_raw_95*tau*low_res_observed_std;
    low_res_lower_contour_raw_95_volume_prct = sum(low_res_lower_contour_raw_95(:))/low_res_middle_contour_volume;
    low_res_upper_contour_raw_95_volume_prct = sum(low_res_upper_contour_raw_95(:))/low_res_middle_contour_volume;
    low_res_mid_on_upper_raw_95              = low_res_upper_contour_raw_95.*low_res_middle_contour;
    low_res_lower_on_mid_raw_95              = low_res_middle_contour.*low_res_lower_contour_raw_95;
    low_res_upper_subset_mid_raw_95          = low_res_upper_contour_raw_95 - low_res_mid_on_upper_raw_95;
    low_res_mid_subset_lower_raw_95          = low_res_middle_contour - low_res_lower_on_mid_raw_95;  

    % Linear Boundary
    low_res_supGa_raw_80_linear                     = prctile(low_res_supG_raw_linear, 80);
    low_res_supGa_raw_90_linear                     = prctile(low_res_supG_raw_linear, 90);
    low_res_supGa_raw_95_linear                     = prctile(low_res_supG_raw_linear, 95);
       
    low_res_lower_contour_raw_80_linear             = low_res_observed_mean >= thr - low_res_supGa_raw_80_linear*tau*low_res_observed_std;
    low_res_upper_contour_raw_80_linear             = low_res_observed_mean >= thr + low_res_supGa_raw_80_linear*tau*low_res_observed_std;
    low_res_lower_contour_raw_80_linear_volume_prct = sum(low_res_lower_contour_raw_80_linear(:))/low_res_middle_contour_volume;
    low_res_upper_contour_raw_80_linear_volume_prct = sum(low_res_upper_contour_raw_80_linear(:))/low_res_middle_contour_volume;
    low_res_mid_on_upper_raw_80_linear              = low_res_upper_contour_raw_80_linear.*low_res_middle_contour;
    low_res_lower_on_mid_raw_80_linear              = low_res_middle_contour.*low_res_lower_contour_raw_80_linear;
    low_res_upper_subset_mid_raw_80_linear          = low_res_upper_contour_raw_80_linear - low_res_mid_on_upper_raw_80_linear;
    low_res_mid_subset_lower_raw_80_linear          = low_res_middle_contour - low_res_lower_on_mid_raw_80_linear;
    
    low_res_lower_contour_raw_90_linear             = low_res_observed_mean >= thr - low_res_supGa_raw_90_linear*tau*low_res_observed_std;
    low_res_upper_contour_raw_90_linear             = low_res_observed_mean >= thr + low_res_supGa_raw_90_linear*tau*low_res_observed_std;
    low_res_lower_contour_raw_90_linear_volume_prct = sum(low_res_lower_contour_raw_90_linear(:))/low_res_middle_contour_volume;
    low_res_upper_contour_raw_90_linear_volume_prct = sum(low_res_upper_contour_raw_90_linear(:))/low_res_middle_contour_volume;
    low_res_mid_on_upper_raw_90_linear              = low_res_upper_contour_raw_90_linear.*low_res_middle_contour;
    low_res_lower_on_mid_raw_90_linear              = low_res_middle_contour.*low_res_lower_contour_raw_90_linear;
    low_res_upper_subset_mid_raw_90_linear          = low_res_upper_contour_raw_90_linear - low_res_mid_on_upper_raw_90_linear;
    low_res_mid_subset_lower_raw_90_linear          = low_res_middle_contour - low_res_lower_on_mid_raw_90_linear;
    
    low_res_lower_contour_raw_95_linear             = low_res_observed_mean >= thr - low_res_supGa_raw_95_linear*tau*low_res_observed_std;
    low_res_upper_contour_raw_95_linear             = low_res_observed_mean >= thr + low_res_supGa_raw_95_linear*tau*low_res_observed_std;
    low_res_lower_contour_raw_95_linear_volume_prct = sum(low_res_lower_contour_raw_95_linear(:))/low_res_middle_contour_volume;
    low_res_upper_contour_raw_95_linear_volume_prct = sum(low_res_upper_contour_raw_95_linear(:))/low_res_middle_contour_volume;
    low_res_mid_on_upper_raw_95_linear              = low_res_upper_contour_raw_95_linear.*low_res_middle_contour;
    low_res_lower_on_mid_raw_95_linear              = low_res_middle_contour.*low_res_lower_contour_raw_95_linear;
    low_res_upper_subset_mid_raw_95_linear          = low_res_upper_contour_raw_95_linear - low_res_mid_on_upper_raw_95_linear;
    low_res_mid_subset_lower_raw_95_linear          = low_res_middle_contour - low_res_lower_on_mid_raw_95_linear;

    %
    % Storing all variables of interest
    %
    low_res_supG_raw_store(:,t)                                    = low_res_supG_raw;
    low_res_threshold_raw_80_store(t)                              = low_res_supGa_raw_80;
    low_res_lower_contour_raw_80_store(t,:,:)                      = low_res_lower_contour_raw_80;
    low_res_upper_contour_raw_80_store(t,:,:)                      = low_res_upper_contour_raw_80;
    low_res_upper_subset_mid_raw_80_store(t,:,:)                   = low_res_upper_subset_mid_raw_80;
    low_res_mid_subset_lower_raw_80_store(t,:,:)                   = low_res_mid_subset_lower_raw_80;
    low_res_lower_contour_raw_80_volume_prct_store(t)              = low_res_lower_contour_raw_80_volume_prct;
    low_res_upper_contour_raw_80_volume_prct_store(t)              = low_res_upper_contour_raw_80_volume_prct;
 
    low_res_threshold_raw_90_store(t)                              = low_res_supGa_raw_90;
    low_res_lower_contour_raw_90_store(t,:,:)                      = low_res_lower_contour_raw_90;
    low_res_upper_contour_raw_90_store(t,:,:)                      = low_res_upper_contour_raw_90;
    low_res_upper_subset_mid_raw_90_store(t,:,:)                   = low_res_upper_subset_mid_raw_90;
    low_res_mid_subset_lower_raw_90_store(t,:,:)                   = low_res_mid_subset_lower_raw_90;
    low_res_lower_contour_raw_90_volume_prct_store(t)              = low_res_lower_contour_raw_90_volume_prct;
    low_res_upper_contour_raw_90_volume_prct_store(t)              = low_res_upper_contour_raw_90_volume_prct;

    low_res_threshold_raw_95_store(t)                              = low_res_supGa_raw_95;
    low_res_lower_contour_raw_95_store(t,:,:)                      = low_res_lower_contour_raw_95;
    low_res_upper_contour_raw_95_store(t,:,:)                      = low_res_upper_contour_raw_95;
    low_res_upper_subset_mid_raw_95_store(t,:,:)                   = low_res_upper_subset_mid_raw_95;
    low_res_mid_subset_lower_raw_95_store(t,:,:)                   = low_res_mid_subset_lower_raw_95;
    low_res_lower_contour_raw_95_volume_prct_store(t)              = low_res_lower_contour_raw_95_volume_prct;
    low_res_upper_contour_raw_95_volume_prct_store(t)              = low_res_upper_contour_raw_95_volume_prct;

    low_res_supG_raw_linear_store(:,t)                                    = low_res_supG_raw_linear;
    low_res_threshold_raw_80_linear_store(t)                              = low_res_supGa_raw_80_linear;
    low_res_lower_contour_raw_80_linear_store(t,:,:)                      = low_res_lower_contour_raw_80_linear;
    low_res_upper_contour_raw_80_linear_store(t,:,:)                      = low_res_upper_contour_raw_80_linear;
    low_res_upper_subset_mid_raw_80_linear_store(t,:,:)                   = low_res_upper_subset_mid_raw_80_linear;
    low_res_mid_subset_lower_raw_80_linear_store(t,:,:)                   = low_res_mid_subset_lower_raw_80_linear;
    low_res_lower_contour_raw_80_linear_volume_prct_store(t)              = low_res_lower_contour_raw_80_linear_volume_prct;
    low_res_upper_contour_raw_80_linear_volume_prct_store(t)              = low_res_upper_contour_raw_80_linear_volume_prct;
 
    low_res_threshold_raw_90_linear_store(t)                              = low_res_supGa_raw_90_linear;
    low_res_lower_contour_raw_90_linear_store(t,:,:)                      = low_res_lower_contour_raw_90_linear;
    low_res_upper_contour_raw_90_linear_store(t,:,:)                      = low_res_upper_contour_raw_90_linear;
    low_res_upper_subset_mid_raw_90_linear_store(t,:,:)                   = low_res_upper_subset_mid_raw_90_linear;
    low_res_mid_subset_lower_raw_90_linear_store(t,:,:)                   = low_res_mid_subset_lower_raw_90_linear;
    low_res_lower_contour_raw_90_linear_volume_prct_store(t)              = low_res_lower_contour_raw_90_linear_volume_prct;
    low_res_upper_contour_raw_90_linear_volume_prct_store(t)              = low_res_upper_contour_raw_90_linear_volume_prct;

    low_res_threshold_raw_95_linear_store(t)                              = low_res_supGa_raw_95_linear;
    low_res_lower_contour_raw_95_linear_store(t,:,:)                      = low_res_lower_contour_raw_95_linear;
    low_res_upper_contour_raw_95_linear_store(t,:,:)                      = low_res_upper_contour_raw_95_linear;
    low_res_upper_subset_mid_raw_95_linear_store(t,:,:)                   = low_res_upper_subset_mid_raw_95_linear;
    low_res_mid_subset_lower_raw_95_linear_store(t,:,:)                   = low_res_mid_subset_lower_raw_95_linear;
    low_res_lower_contour_raw_95_linear_volume_prct_store(t)              = low_res_lower_contour_raw_95_linear_volume_prct;
    low_res_upper_contour_raw_95_linear_volume_prct_store(t)              = low_res_upper_contour_raw_95_linear_volume_prct;
    
    
    if sum(low_res_upper_subset_mid_raw_80(:))+sum(low_res_mid_subset_lower_raw_80(:))==0
      low_res_subset_success_vector_raw_80(t) = 1; 
      fprintf('low res raw nominal 80 true boundary success! \n');
    else 
      low_res_subset_success_vector_raw_80(t) = 0; 
      fprintf('low res raw nominal 80 true boundary failure! \n');
    end 

    if sum(low_res_upper_subset_mid_raw_90(:))+sum(low_res_mid_subset_lower_raw_90(:))==0
      low_res_subset_success_vector_raw_90(t) = 1; 
      fprintf('low res raw nominal 90 true boundary success! \n');
    else 
      low_res_subset_success_vector_raw_90(t) = 0; 
      fprintf('low res raw nominal 90 true boundary failure! \n');
    end 

    if sum(low_res_upper_subset_mid_raw_95(:))+sum(low_res_mid_subset_lower_raw_95(:))==0
      low_res_subset_success_vector_raw_95(t) = 1; 
      fprintf('low res raw nominal 95 true boundary success! \n');
    else 
      low_res_subset_success_vector_raw_95(t) = 0; 
      fprintf('low res raw nominal 95 true boundary failure! \n');
    end 

    if sum(low_res_upper_subset_mid_raw_80_linear(:))+sum(low_res_mid_subset_lower_raw_80_linear(:))==0
      low_res_subset_success_vector_raw_80_linear(t) = 1; 
      fprintf('low res raw nominal 80 linear boundary success! \n');
    else 
      low_res_subset_success_vector_raw_80_linear(t) = 0; 
      fprintf('low res raw nominal 80 linear boundary failure! \n');
    end 

    if sum(low_res_upper_subset_mid_raw_90_linear(:))+sum(low_res_mid_subset_lower_raw_90_linear(:))==0
      low_res_subset_success_vector_raw_90_linear(t) = 1; 
      fprintf('low res raw nominal 90 linear boundary success! \n');
    else 
      low_res_subset_success_vector_raw_90_linear(t) = 0; 
      fprintf('low res raw nominal 90 linear boundary failure! \n');
    end 

    if sum(low_res_upper_subset_mid_raw_95_linear(:))+sum(low_res_mid_subset_lower_raw_95_linear(:))==0
      low_res_subset_success_vector_raw_95_linear(t) = 1; 
      fprintf('low res raw nominal 95 linear boundary success! \n');
    else 
      low_res_subset_success_vector_raw_95_linear(t) = 0; 
      fprintf('low res raw nominal 95 linear boundary failure! \n');
    end
    
end

percentage_success_vector_raw_80                         = mean(subset_success_vector_raw_80, 1);
percentage_success_vector_raw_90                         = mean(subset_success_vector_raw_90, 1);
percentage_success_vector_raw_95                         = mean(subset_success_vector_raw_95, 1);

percentage_success_vector_raw_80_linear                  = mean(subset_success_vector_raw_80_linear, 1);
percentage_success_vector_raw_90_linear                  = mean(subset_success_vector_raw_90_linear, 1);
percentage_success_vector_raw_95_linear                  = mean(subset_success_vector_raw_95_linear, 1);

low_res_percentage_success_vector_raw_80                         = mean(low_res_subset_success_vector_raw_80, 1);
low_res_percentage_success_vector_raw_90                         = mean(low_res_subset_success_vector_raw_90, 1);
low_res_percentage_success_vector_raw_95                         = mean(low_res_subset_success_vector_raw_95, 1);

low_res_percentage_success_vector_raw_80_linear                  = mean(low_res_subset_success_vector_raw_80_linear, 1);
low_res_percentage_success_vector_raw_90_linear                  = mean(low_res_subset_success_vector_raw_90_linear, 1);
low_res_percentage_success_vector_raw_95_linear                  = mean(low_res_subset_success_vector_raw_95_linear, 1);


eval(['save ' SvNm ' nSubj nRlz dim smo mag rimFWHM thr nBoot '... 
      'threshold_raw_80_store threshold_raw_90_store threshold_raw_95_store threshold_raw_80_linear_store threshold_raw_90_linear_store threshold_raw_95_linear_store low_res_threshold_raw_80_store low_res_threshold_raw_90_store low_res_threshold_raw_95_store low_res_threshold_raw_80_linear_store low_res_threshold_raw_90_linear_store low_res_threshold_raw_95_linear_store '...
      'lower_contour_raw_80_store lower_contour_raw_90_store lower_contour_raw_95_store lower_contour_raw_80_linear_store lower_contour_raw_90_linear_store lower_contour_raw_95_linear_store low_res_lower_contour_raw_80_store low_res_lower_contour_raw_90_store low_res_lower_contour_raw_95_store low_res_lower_contour_raw_80_linear_store low_res_lower_contour_raw_90_linear_store low_res_lower_contour_raw_95_linear_store '...
      'upper_contour_raw_80_store upper_contour_raw_90_store upper_contour_raw_95_store upper_contour_raw_80_linear_store upper_contour_raw_90_linear_store upper_contour_raw_95_linear_store low_res_upper_contour_raw_80_store low_res_upper_contour_raw_90_store low_res_upper_contour_raw_95_store low_res_upper_contour_raw_80_linear_store low_res_upper_contour_raw_90_linear_store low_res_upper_contour_raw_95_linear_store '...
      'upper_subset_mid_raw_80_store upper_subset_mid_raw_90_store upper_subset_mid_raw_95_store upper_subset_mid_raw_80_linear_store upper_subset_mid_raw_90_linear_store upper_subset_mid_raw_95_linear_store low_res_upper_subset_mid_raw_80_store low_res_upper_subset_mid_raw_90_store low_res_upper_subset_mid_raw_95_store low_res_upper_subset_mid_raw_80_linear_store low_res_upper_subset_mid_raw_90_linear_store low_res_upper_subset_mid_raw_95_linear_store '...
      'mid_subset_lower_raw_80_store mid_subset_lower_raw_90_store mid_subset_lower_raw_95_store mid_subset_lower_raw_80_linear_store mid_subset_lower_raw_90_linear_store mid_subset_lower_raw_95_linear_store low_res_mid_subset_lower_raw_80_store low_res_mid_subset_lower_raw_90_store low_res_mid_subset_lower_raw_95_store low_res_mid_subset_lower_raw_80_linear_store low_res_mid_subset_lower_raw_90_linear_store low_res_mid_subset_lower_raw_95_linear_store '...
      'subset_success_vector_raw_80 subset_success_vector_raw_90 subset_success_vector_raw_95 subset_success_vector_raw_80_linear subset_success_vector_raw_90_linear subset_success_vector_raw_95_linear low_res_subset_success_vector_raw_80 low_res_subset_success_vector_raw_90 low_res_subset_success_vector_raw_95 low_res_subset_success_vector_raw_80_linear low_res_subset_success_vector_raw_90_linear low_res_subset_success_vector_raw_95_linear '...
      'percentage_success_vector_raw_80 percentage_success_vector_raw_90 percentage_success_vector_raw_95 percentage_success_vector_raw_80_linear percentage_success_vector_raw_90_linear percentage_success_vector_raw_95_linear low_res_percentage_success_vector_raw_80 low_res_percentage_success_vector_raw_90 low_res_percentage_success_vector_raw_95 low_res_percentage_success_vector_raw_80_linear low_res_percentage_success_vector_raw_90_linear low_res_percentage_success_vector_raw_95_linear '...
      'supG_raw_store supG_raw_linear_store low_res_supG_raw_store low_res_supG_raw_linear_store '...
      'middle_contour_volume observed_AC_volume low_res_middle_contour_volume low_res_observed_AC_volume '...
      'lower_contour_raw_80_volume_prct_store lower_contour_raw_90_volume_prct_store lower_contour_raw_95_volume_prct_store lower_contour_raw_80_linear_volume_prct_store lower_contour_raw_90_linear_volume_prct_store lower_contour_raw_95_linear_volume_prct_store low_res_lower_contour_raw_80_volume_prct_store low_res_lower_contour_raw_90_volume_prct_store low_res_lower_contour_raw_95_volume_prct_store low_res_lower_contour_raw_80_linear_volume_prct_store low_res_lower_contour_raw_90_linear_volume_prct_store low_res_lower_contour_raw_95_linear_volume_prct_store '...
      'upper_contour_raw_80_volume_prct_store upper_contour_raw_90_volume_prct_store upper_contour_raw_95_volume_prct_store upper_contour_raw_80_linear_volume_prct_store upper_contour_raw_90_linear_volume_prct_store upper_contour_raw_95_linear_volume_prct_store low_res_upper_contour_raw_80_volume_prct_store low_res_upper_contour_raw_90_volume_prct_store low_res_upper_contour_raw_95_volume_prct_store low_res_upper_contour_raw_80_linear_volume_prct_store low_res_upper_contour_raw_90_linear_volume_prct_store low_res_upper_contour_raw_95_linear_volume_prct_store'])
