function Sim_15(nSubj,SvNm,nRlz)
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
dim     = [100 100 100]; 
mag     = 3;
rad     = 30;
smo     = 10;
rimFWHM = 15;
stdblk  = prod(dim([1 2])/2);
thr     = 2;

%-----------Initialization of Some Variables
V           = prod(dim);   
wdim        = dim + 2*ceil(rimFWHM*smo)*ones(1,3);  % Working image dimension
trunc_x     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(1))};
trunc_y     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(2))};
trunc_z     = {(ceil(rimFWHM*smo)+1):(ceil(rimFWHM*smo)+dim(3))};
trnind      = cat(2, trunc_x, trunc_y, trunc_z);

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

%- This vector stores the threshold value 'c' for each run
threshold_raw_80_store                  = zeros(nRlz, 1);
threshold_raw_90_store                  = zeros(nRlz, 1);
threshold_raw_95_store                  = zeros(nRlz, 1);

threshold_raw_80_linear_store           = zeros(nRlz, 1);
threshold_raw_90_linear_store           = zeros(nRlz, 1);
threshold_raw_95_linear_store           = zeros(nRlz, 1);

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

observed_middle_contour_volume_store            = zeros(nRlz, 1);

% This stores the vector SupG for each run
supG_raw_store         = zeros(nBoot, nRlz);
supG_raw_linear_store  = zeros(nBoot, nRlz);

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

lower_contour_raw_90_store                       = zeros([nRlz dim]);
upper_contour_raw_90_store                       = zeros([nRlz dim]);
upper_subset_mid_raw_90_store                    = zeros([nRlz dim]);
mid_subset_lower_raw_90_store                    = zeros([nRlz dim]);
lower_contour_raw_90_linear_store                = zeros([nRlz dim]);
upper_contour_raw_90_linear_store                = zeros([nRlz dim]);
upper_subset_mid_raw_90_linear_store             = zeros([nRlz dim]);
mid_subset_lower_raw_90_linear_store             = zeros([nRlz dim]);

lower_contour_raw_95_store                       = zeros([nRlz dim]);
upper_contour_raw_95_store                       = zeros([nRlz dim]);
upper_subset_mid_raw_95_store                    = zeros([nRlz dim]);
mid_subset_lower_raw_95_store                    = zeros([nRlz dim]);
lower_contour_raw_95_linear_store                = zeros([nRlz dim]);
upper_contour_raw_95_linear_store                = zeros([nRlz dim]);
upper_subset_mid_raw_95_linear_store             = zeros([nRlz dim]);
mid_subset_lower_raw_95_linear_store             = zeros([nRlz dim]);

supG_raw              = zeros(nBoot,1);
supG_raw_linear       = zeros(nBoot,1);

% Creating linearly increasing signal across columns
Sig = SpheroidSignal(wdim, rad, mag, 0);

% Smoothing the signal
Sigs = zeros(wdim);
ss   = spm_smooth(Sig,Sigs,smo*ones(1,3));
Sigs = Sigs;

% Truncate to avoid edge effects
tSigs          = Sigs(trnind{1}, trnind{2}, trnind{3});
maxtSigs       = max(tSigs(:));
Sig            = (mag/maxtSigs)*tSigs;

% Uncomment to look at the Signal
%imagesc(Sig); axis image; colorbar
AC = Sig >= thr;
middle_contour                = AC;
middle_contour_volume         = sum(middle_contour(:));

% Variables for computing the estimated boundary
[a,b,c] = ndgrid(-1:1);
se = strel('arbitrary',sqrt(a.^2 + b.^2 + c.^2) <=1);

for t=1:nRlz
    fprintf('.');
      for i=1:nSubj
	    %
	    % Generate random realizations of signal + noise
	    %
        raw_noise(:,:,:,i) = randn(wdim); %- Noise that will be added to the signal 

        %
        % smooth noise  
        %
        Noises = zeros(wdim);
        tt     = spm_smooth(raw_noise(:,:,:,i),Noises,smo*ones(1,3));
        Noises = Noises/sqrt(tt);      
      
        %
        % Truncate to avoid edge effects
        %
        tNoises = Noises(trnind{1},trnind{2},trnind{3});       
        tImgs = Sig + tNoises; % Creates the true image of smoothed signal + smoothed noise
        observed_data(:,:,:,i) = tImgs;
        
      end %========== Loop i (subjects)
      
      observed_mean = mean(observed_data,4);

      observed_std = reshape(...
         biasmystd(reshape(observed_data,[prod(dim) nSubj]),stdblk),...
           dim);
       
      % Making the three observed boundaries: dilated boundary, eroded
      % boundary, and dilated - eroded boundary.
      observed_AC = observed_mean >= thr;
      observed_AC_volume = sum(observed_AC(:)); 
 
      % Residuals
      resid = bsxfun(@minus,observed_data,observed_mean);
      resid = spdiags(1./reshape(observed_std, [prod(dim) 1]), 0,prod(dim),prod(dim))*reshape(resid,[prod(dim) nSubj]);
      resid = reshape(resid, [dim nSubj]);
      
      % Extracting only the residuals along the true boundary (using the
      % linear interpolation method). 
      resid_true_boundary = linear_interp_boundary(AC, resid, nSubj);
      
      % Extracting only the residuals along the esitimated boundary (using the linear interpolation method)
      resid_estimated_boundary = linear_interp_boundary(observed_AC, resid, nSubj);
      
      %% Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          % Applying the bootstrap using Rademacher variables (signflips)
          signflips                              = randi(2,[nSubj,1])*2-3;
          
          % True boundary
          true_boundary_bootstrap                = resid_true_boundary*spdiags(signflips, 0, nSubj, nSubj);
          true_boundary_resid_field              = sum(true_boundary_bootstrap, 2)/sqrt(nSubj); 
          supG_raw(k)                            = max(abs(true_boundary_resid_field));
          
          % Estimated boundary
          estimated_boundary_bootstrap           = resid_estimated_boundary*spdiags(signflips, 0, nSubj, nSubj);
          estimated_boundary_resid_field         = sum(estimated_boundary_bootstrap, 2)/sqrt(nSubj); 
          supG_raw_linear(k)                     = max(abs(estimated_boundary_resid_field));
          
      end
    
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
    supG_raw_store(:,t)                                      = supG_raw;
    threshold_raw_80_store(t)                                = supGa_raw_80;
    lower_contour_raw_80_store(t,:,:,:)                      = lower_contour_raw_80;
    upper_contour_raw_80_store(t,:,:,:)                      = upper_contour_raw_80;
    upper_subset_mid_raw_80_store(t,:,:,:)                   = upper_subset_mid_raw_80;
    mid_subset_lower_raw_80_store(t,:,:,:)                   = mid_subset_lower_raw_80;
    lower_contour_raw_80_volume_prct_store(t)                = lower_contour_raw_80_volume_prct;
    upper_contour_raw_80_volume_prct_store(t)                = upper_contour_raw_80_volume_prct;
 
    threshold_raw_90_store(t)                                = supGa_raw_90;
    lower_contour_raw_90_store(t,:,:,:)                      = lower_contour_raw_90;
    upper_contour_raw_90_store(t,:,:,:)                      = upper_contour_raw_90;
    upper_subset_mid_raw_90_store(t,:,:,:)                   = upper_subset_mid_raw_90;
    mid_subset_lower_raw_90_store(t,:,:,:)                   = mid_subset_lower_raw_90;
    lower_contour_raw_90_volume_prct_store(t)                = lower_contour_raw_90_volume_prct;
    upper_contour_raw_90_volume_prct_store(t)                = upper_contour_raw_90_volume_prct;

    threshold_raw_95_store(t)                                = supGa_raw_95;
    lower_contour_raw_95_store(t,:,:,:)                      = lower_contour_raw_95;
    upper_contour_raw_95_store(t,:,:,:)                      = upper_contour_raw_95;
    upper_subset_mid_raw_95_store(t,:,:,:)                   = upper_subset_mid_raw_95;
    mid_subset_lower_raw_95_store(t,:,:,:)                   = mid_subset_lower_raw_95;
    lower_contour_raw_95_volume_prct_store(t)                = lower_contour_raw_95_volume_prct;
    upper_contour_raw_95_volume_prct_store(t)                = upper_contour_raw_95_volume_prct;

    supG_raw_linear_store(:,t)                                      = supG_raw_linear;
    threshold_raw_80_linear_store(t)                                = supGa_raw_80_linear;
    lower_contour_raw_80_linear_store(t,:,:,:)                      = lower_contour_raw_80_linear;
    upper_contour_raw_80_linear_store(t,:,:,:)                      = upper_contour_raw_80_linear;
    upper_subset_mid_raw_80_linear_store(t,:,:,:)                   = upper_subset_mid_raw_80_linear;
    mid_subset_lower_raw_80_linear_store(t,:,:,:)                   = mid_subset_lower_raw_80_linear;
    lower_contour_raw_80_linear_volume_prct_store(t)                = lower_contour_raw_80_linear_volume_prct;
    upper_contour_raw_80_linear_volume_prct_store(t)                = upper_contour_raw_80_linear_volume_prct;
 
    threshold_raw_90_linear_store(t)                                = supGa_raw_90_linear;
    lower_contour_raw_90_linear_store(t,:,:,:)                      = lower_contour_raw_90_linear;
    upper_contour_raw_90_linear_store(t,:,:,:)                      = upper_contour_raw_90_linear;
    upper_subset_mid_raw_90_linear_store(t,:,:,:)                   = upper_subset_mid_raw_90_linear;
    mid_subset_lower_raw_90_linear_store(t,:,:,:)                   = mid_subset_lower_raw_90_linear;
    lower_contour_raw_90_linear_volume_prct_store(t)                = lower_contour_raw_90_linear_volume_prct;
    upper_contour_raw_90_linear_volume_prct_store(t)                = upper_contour_raw_90_linear_volume_prct;

    threshold_raw_95_linear_store(t)                                = supGa_raw_95_linear;
    lower_contour_raw_95_linear_store(t,:,:,:)                      = lower_contour_raw_95_linear;
    upper_contour_raw_95_linear_store(t,:,:,:)                      = upper_contour_raw_95_linear;
    upper_subset_mid_raw_95_linear_store(t,:,:,:)                   = upper_subset_mid_raw_95_linear;
    mid_subset_lower_raw_95_linear_store(t,:,:,:)                   = mid_subset_lower_raw_95_linear;
    lower_contour_raw_95_linear_volume_prct_store(t)                = lower_contour_raw_95_linear_volume_prct;
    upper_contour_raw_95_linear_volume_prct_store(t)                = upper_contour_raw_95_linear_volume_prct;
    
    observed_middle_contour_volume_store(t)                         = observed_AC_volume;
    
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

end

percentage_success_vector_raw_80                         = mean(subset_success_vector_raw_80, 1);
percentage_success_vector_raw_90                         = mean(subset_success_vector_raw_90, 1);
percentage_success_vector_raw_95                         = mean(subset_success_vector_raw_95, 1);

percentage_success_vector_raw_80_linear                  = mean(subset_success_vector_raw_80_linear, 1);
percentage_success_vector_raw_90_linear                  = mean(subset_success_vector_raw_90_linear, 1);
percentage_success_vector_raw_95_linear                  = mean(subset_success_vector_raw_95_linear, 1);

eval(['save ' SvNm ' nSubj nRlz dim smo mag rimFWHM thr nBoot '... 
      'threshold_raw_80_store threshold_raw_90_store threshold_raw_95_store threshold_raw_80_linear_store threshold_raw_90_linear_store threshold_raw_95_linear_store low_res_threshold_raw_80_store low_res_threshold_raw_90_store low_res_threshold_raw_95_store low_res_threshold_raw_80_linear_store low_res_threshold_raw_90_linear_store low_res_threshold_raw_95_linear_store '...
      'lower_contour_raw_80_store lower_contour_raw_90_store lower_contour_raw_95_store lower_contour_raw_80_linear_store lower_contour_raw_90_linear_store lower_contour_raw_95_linear_store low_res_lower_contour_raw_80_store low_res_lower_contour_raw_90_store low_res_lower_contour_raw_95_store low_res_lower_contour_raw_80_linear_store low_res_lower_contour_raw_90_linear_store low_res_lower_contour_raw_95_linear_store '...
      'upper_contour_raw_80_store upper_contour_raw_90_store upper_contour_raw_95_store upper_contour_raw_80_linear_store upper_contour_raw_90_linear_store upper_contour_raw_95_linear_store low_res_upper_contour_raw_80_store low_res_upper_contour_raw_90_store low_res_upper_contour_raw_95_store low_res_upper_contour_raw_80_linear_store low_res_upper_contour_raw_90_linear_store low_res_upper_contour_raw_95_linear_store '...
      'upper_subset_mid_raw_80_store upper_subset_mid_raw_90_store upper_subset_mid_raw_95_store upper_subset_mid_raw_80_linear_store upper_subset_mid_raw_90_linear_store upper_subset_mid_raw_95_linear_store low_res_upper_subset_mid_raw_80_store low_res_upper_subset_mid_raw_90_store low_res_upper_subset_mid_raw_95_store low_res_upper_subset_mid_raw_80_linear_store low_res_upper_subset_mid_raw_90_linear_store low_res_upper_subset_mid_raw_95_linear_store '...
      'mid_subset_lower_raw_80_store mid_subset_lower_raw_90_store mid_subset_lower_raw_95_store mid_subset_lower_raw_80_linear_store mid_subset_lower_raw_90_linear_store mid_subset_lower_raw_95_linear_store low_res_mid_subset_lower_raw_80_store low_res_mid_subset_lower_raw_90_store low_res_mid_subset_lower_raw_95_store low_res_mid_subset_lower_raw_80_linear_store low_res_mid_subset_lower_raw_90_linear_store low_res_mid_subset_lower_raw_95_linear_store '...
      'subset_success_vector_raw_80 subset_success_vector_raw_90 subset_success_vector_raw_95 subset_success_vector_raw_80_linear subset_success_vector_raw_90_linear subset_success_vector_raw_95_linear low_res_subset_success_vector_raw_80 low_res_subset_success_vector_raw_90 low_res_subset_success_vector_raw_95 low_res_subset_success_vector_raw_80_linear low_res_subset_success_vector_raw_90_linear low_res_subset_success_vector_raw_95_linear '...
      'percentage_success_vector_raw_80 percentage_success_vector_raw_90 percentage_success_vector_raw_95 percentage_success_vector_raw_80_linear percentage_success_vector_raw_90_linear percentage_success_vector_raw_95_linear low_res_percentage_success_vector_raw_80 low_res_percentage_success_vector_raw_90 low_res_percentage_success_vector_raw_95 low_res_percentage_success_vector_raw_80_linear low_res_percentage_success_vector_raw_90_linear low_res_percentage_success_vector_raw_95_linear '...
      'supG_raw_store supG_raw_linear_store low_res_supG_raw_store low_res_supG_raw_linear_store '...
      'middle_contour middle_contour_volume observed_middle_contour_volume low_res_middle_contour_volume low_res_observed_AC_volume '...
      'lower_contour_raw_80_volume_prct_store lower_contour_raw_90_volume_prct_store lower_contour_raw_95_volume_prct_store lower_contour_raw_80_linear_volume_prct_store lower_contour_raw_90_linear_volume_prct_store lower_contour_raw_95_linear_volume_prct_store low_res_lower_contour_raw_80_volume_prct_store low_res_lower_contour_raw_90_volume_prct_store low_res_lower_contour_raw_95_volume_prct_store low_res_lower_contour_raw_80_linear_volume_prct_store low_res_lower_contour_raw_90_linear_volume_prct_store low_res_lower_contour_raw_95_linear_volume_prct_store '...
      'upper_contour_raw_80_volume_prct_store upper_contour_raw_90_volume_prct_store upper_contour_raw_95_volume_prct_store upper_contour_raw_80_linear_volume_prct_store upper_contour_raw_90_linear_volume_prct_store upper_contour_raw_95_linear_volume_prct_store low_res_upper_contour_raw_80_volume_prct_store low_res_upper_contour_raw_90_volume_prct_store low_res_upper_contour_raw_95_volume_prct_store low_res_upper_contour_raw_80_linear_volume_prct_store low_res_upper_contour_raw_90_linear_volume_prct_store low_res_upper_contour_raw_95_linear_volume_prct_store'])
