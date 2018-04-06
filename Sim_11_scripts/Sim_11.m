function Sim_11(nSubj,SvNm,nRlz)
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
subset_success_vector_raw_80_gaussian           = zeros(nRlz, 1); 
subset_success_vector_raw_90_gaussian           = zeros(nRlz, 1);
subset_success_vector_raw_95_gaussian           = zeros(nRlz, 1); 
subset_success_vector_raw_80_ero_dil_gaussian   = zeros(nRlz, 1); 
subset_success_vector_raw_90_ero_dil_gaussian   = zeros(nRlz, 1);
subset_success_vector_raw_95_ero_dil_gaussian   = zeros(nRlz, 1); 

subset_success_vector_raw_80_signflips           = zeros(nRlz, 1); 
subset_success_vector_raw_90_signflips           = zeros(nRlz, 1);
subset_success_vector_raw_95_signflips           = zeros(nRlz, 1); 
subset_success_vector_raw_80_ero_dil_signflips   = zeros(nRlz, 1); 
subset_success_vector_raw_90_ero_dil_signflips   = zeros(nRlz, 1);
subset_success_vector_raw_95_ero_dil_signflips   = zeros(nRlz, 1); 

%- This vector stores the threshold value 'c' for each run
threshold_raw_80_gaussian_store                  = zeros(nRlz, 1);
threshold_raw_90_gaussian_store                  = zeros(nRlz, 1);
threshold_raw_95_gaussian_store                  = zeros(nRlz, 1);
threshold_raw_80_ero_dil_gaussian_store          = zeros(nRlz, 1);
threshold_raw_90_ero_dil_gaussian_store          = zeros(nRlz, 1);
threshold_raw_95_ero_dil_gaussian_store          = zeros(nRlz, 1);

threshold_raw_80_signflips_store                  = zeros(nRlz, 1);
threshold_raw_90_signflips_store                  = zeros(nRlz, 1);
threshold_raw_95_signflips_store                  = zeros(nRlz, 1);
threshold_raw_80_ero_dil_signflips_store          = zeros(nRlz, 1);
threshold_raw_90_ero_dil_signflips_store          = zeros(nRlz, 1);
threshold_raw_95_ero_dil_signflips_store          = zeros(nRlz, 1);

%- This vector stores the percentage volumes A^+_c/A_c, A^_c/A_c, A^-_c/A_c
lower_contour_raw_80_gaussian_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_80_gaussian_volume_prct_store = zeros(nRlz, 1);
lower_contour_raw_80_ero_dil_gaussian_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_80_ero_dil_gaussian_volume_prct_store = zeros(nRlz, 1);
lower_contour_raw_80_signflips_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_80_signflips_volume_prct_store = zeros(nRlz, 1);
lower_contour_raw_80_ero_dil_signflips_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_80_ero_dil_signflips_volume_prct_store = zeros(nRlz, 1);

lower_contour_raw_90_gaussian_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_90_gaussian_volume_prct_store = zeros(nRlz, 1);
lower_contour_raw_90_ero_dil_gaussian_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_90_ero_dil_gaussian_volume_prct_store = zeros(nRlz, 1);
lower_contour_raw_90_signflips_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_90_signflips_volume_prct_store = zeros(nRlz, 1);
lower_contour_raw_90_ero_dil_signflips_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_90_ero_dil_signflips_volume_prct_store = zeros(nRlz, 1);

lower_contour_raw_95_gaussian_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_95_gaussian_volume_prct_store = zeros(nRlz, 1);
lower_contour_raw_95_ero_dil_gaussian_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_95_ero_dil_gaussian_volume_prct_store = zeros(nRlz, 1);
lower_contour_raw_95_signflips_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_95_signflips_volume_prct_store = zeros(nRlz, 1);
lower_contour_raw_95_ero_dil_signflips_volume_prct_store = zeros(nRlz, 1);
upper_contour_raw_95_ero_dil_signflips_volume_prct_store = zeros(nRlz, 1);

% This stores the vector SupG for each run
supG_raw_gaussian_store          = zeros(nBoot, nRlz);
supG_raw_ero_dil_gaussian_store  = zeros(nBoot, nRlz);
supG_raw_signflips_store         = zeros(nBoot, nRlz);
supG_raw_ero_dil_signflips_store = zeros(nBoot, nRlz);

%-These matrices store all the sets of interest during the bootstrap
% method for all levels of smoothing
lower_contour_raw_80_gaussian_store               = zeros([nRlz dim]);
upper_contour_raw_80_gaussian_store               = zeros([nRlz dim]);
upper_subset_mid_raw_80_gaussian_store            = zeros([nRlz dim]);
mid_subset_lower_raw_80_gaussian_store            = zeros([nRlz dim]);
lower_contour_raw_80_ero_dil_gaussian_store       = zeros([nRlz dim]);
upper_contour_raw_80_ero_dil_gaussian_store       = zeros([nRlz dim]);
upper_subset_mid_raw_80_ero_dil_gaussian_store    = zeros([nRlz dim]);
mid_subset_lower_raw_80_ero_dil_gaussian_store    = zeros([nRlz dim]);

lower_contour_raw_80_signflips_store               = zeros([nRlz dim]);
upper_contour_raw_80_signflips_store               = zeros([nRlz dim]);
upper_subset_mid_raw_80_signflips_store            = zeros([nRlz dim]);
mid_subset_lower_raw_80_signflips_store            = zeros([nRlz dim]);
lower_contour_raw_80_ero_dil_signflips_store       = zeros([nRlz dim]);
upper_contour_raw_80_ero_dil_signflips_store       = zeros([nRlz dim]);
upper_subset_mid_raw_80_ero_dil_signflips_store    = zeros([nRlz dim]);
mid_subset_lower_raw_80_ero_dil_signflips_store    = zeros([nRlz dim]);

lower_contour_raw_90_gaussian_store               = zeros([nRlz dim]);
upper_contour_raw_90_gaussian_store               = zeros([nRlz dim]);
upper_subset_mid_raw_90_gaussian_store            = zeros([nRlz dim]);
mid_subset_lower_raw_90_gaussian_store            = zeros([nRlz dim]);
lower_contour_raw_90_ero_dil_gaussian_store       = zeros([nRlz dim]);
upper_contour_raw_90_ero_dil_gaussian_store       = zeros([nRlz dim]);
upper_subset_mid_raw_90_ero_dil_gaussian_store    = zeros([nRlz dim]);
mid_subset_lower_raw_90_ero_dil_gaussian_store    = zeros([nRlz dim]);

lower_contour_raw_90_signflips_store               = zeros([nRlz dim]);
upper_contour_raw_90_signflips_store               = zeros([nRlz dim]);
upper_subset_mid_raw_90_signflips_store            = zeros([nRlz dim]);
mid_subset_lower_raw_90_signflips_store            = zeros([nRlz dim]);
lower_contour_raw_90_ero_dil_signflips_store       = zeros([nRlz dim]);
upper_contour_raw_90_ero_dil_signflips_store       = zeros([nRlz dim]);
upper_subset_mid_raw_90_ero_dil_signflips_store    = zeros([nRlz dim]);
mid_subset_lower_raw_90_ero_dil_signflips_store    = zeros([nRlz dim]);

lower_contour_raw_95_gaussian_store               = zeros([nRlz dim]);
upper_contour_raw_95_gaussian_store               = zeros([nRlz dim]);
upper_subset_mid_raw_95_gaussian_store            = zeros([nRlz dim]);
mid_subset_lower_raw_95_gaussian_store            = zeros([nRlz dim]);
lower_contour_raw_95_ero_dil_gaussian_store       = zeros([nRlz dim]);
upper_contour_raw_95_ero_dil_gaussian_store       = zeros([nRlz dim]);
upper_subset_mid_raw_95_ero_dil_gaussian_store    = zeros([nRlz dim]);
mid_subset_lower_raw_95_ero_dil_gaussian_store    = zeros([nRlz dim]);

lower_contour_raw_95_signflips_store               = zeros([nRlz dim]);
upper_contour_raw_95_signflips_store               = zeros([nRlz dim]);
upper_subset_mid_raw_95_signflips_store            = zeros([nRlz dim]);
mid_subset_lower_raw_95_signflips_store            = zeros([nRlz dim]);
lower_contour_raw_95_ero_dil_signflips_store       = zeros([nRlz dim]);
upper_contour_raw_95_ero_dil_signflips_store       = zeros([nRlz dim]);
upper_subset_mid_raw_95_ero_dil_signflips_store    = zeros([nRlz dim]);
mid_subset_lower_raw_95_ero_dil_signflips_store    = zeros([nRlz dim]);

supG_raw_gaussian          = zeros(nBoot,1);
supG_raw_ero_dil_gaussian  = zeros(nBoot,1);
supG_raw_signflips         = zeros(nBoot,1);
supG_raw_ero_dil_signflips = zeros(nBoot,1);


raw_field_boundary_store        = zeros(dim(2), nBoot);

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

      observed_mean = mean(observed_data,3);

      observed_std = reshape(...
         biasmystd(reshape(observed_data,[prod(dim) nSubj]),stdblk),...
           dim);
       
      % Making the three observed boundaries: dilated boundary, eroded
      % boundary, and dilated - eroded boundary.
      observed_AC = observed_mean >= thr;
      observed_AC_volume = sum(observed_AC(:)); 
      observed_AC_ero = imerode(observed_AC,se);
      observed_AC_dil = imdilate(observed_AC,se);
      observed_delta_AC_ero_dil = (observed_AC_dil - observed_AC)|(observed_AC - observed_AC_ero);
       
      % Residuals
      resid = bsxfun(@minus,observed_data,observed_mean);
      resid = spdiags(1./reshape(observed_std, [prod(dim) 1]), 0,prod(dim),prod(dim))*reshape(resid,[prod(dim) nSubj]); 
            
      % Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          % Applying the bootstrap using gaussian random variables (gaussian_rand) and Rademacher variables (signflips)
          gaussian_rand                          = normrnd(0,1,[nSubj, 1]);
          resid_bootstrap_gaussian               = resid*spdiags(gaussian_rand, 0, nSubj, nSubj);
          resid_bootstrap_gaussian               = reshape(resid_bootstrap_gaussian, [dim nSubj]);
          resid_field_gaussian                   = sum(resid_bootstrap_gaussian, 3)/sqrt(nSubj); 

          supG_raw_gaussian(k)          = max(abs(resid_field_gaussian(true_boundary)));
          supG_raw_ero_dil_gaussian(k)  = max(abs(resid_field_gaussian(observed_delta_AC_ero_dil)));

          signflips                               = randi(2,[nSubj,1])*2-3;
          resid_bootstrap_signflips               = resid*spdiags(signflips, 0, nSubj, nSubj);
          resid_bootstrap_signflips               = reshape(resid_bootstrap_signflips, [dim nSubj]);
          resid_field_signflips                   = sum(resid_bootstrap_signflips, 3)/sqrt(nSubj); 

          supG_raw_signflips(k)          = max(abs(resid_field_signflips(true_boundary)));
          supG_raw_ero_dil_signflips(k)  = max(abs(resid_field_signflips(observed_delta_AC_ero_dil)));
          
          if t==nRlz
              raw_field_boundary_store(:,k)         = resid_field_gaussian(true_boundary);
          end 
      end
    middle_contour                = AC;
    middle_contour_volume         = sum(middle_contour(:));
    
    %% Gaussian random variable results for the true and estimated boundary
    % True boundary
    supGa_raw_80_gaussian                     = prctile(supG_raw_gaussian, 80);
    supGa_raw_90_gaussian                     = prctile(supG_raw_gaussian, 90);
    supGa_raw_95_gaussian                     = prctile(supG_raw_gaussian, 95);
       
    lower_contour_raw_80_gaussian             = observed_mean >= thr - supGa_raw_80_gaussian*tau*observed_std;
    upper_contour_raw_80_gaussian             = observed_mean >= thr + supGa_raw_80_gaussian*tau*observed_std;
    lower_contour_raw_80_gaussian_volume_prct = sum(lower_contour_raw_80_gaussian(:))/middle_contour_volume;
    upper_contour_raw_80_gaussian_volume_prct = sum(upper_contour_raw_80_gaussian(:))/middle_contour_volume;
    mid_on_upper_raw_80_gaussian              = upper_contour_raw_80.*middle_contour;
    lower_on_mid_raw_80_gaussian              = middle_contour.*lower_contour_raw_80_gaussian;
    upper_subset_mid_raw_80_gaussian          = upper_contour_raw_80_gaussian - mid_on_upper_raw_80_gaussian;
    mid_subset_lower_raw_80_gaussian          = middle_contour - lower_on_mid_raw_80_gaussian;
    
    lower_contour_raw_90_gaussian             = observed_mean >= thr - supGa_raw_90_gaussian*tau*observed_std;
    upper_contour_raw_90_gaussian             = observed_mean >= thr + supGa_raw_90_gaussian*tau*observed_std;
    lower_contour_raw_90_gaussian_volume_prct = sum(lower_contour_raw_90_gaussian(:))/middle_contour_volume;
    upper_contour_raw_90_gaussian_volume_prct = sum(upper_contour_raw_90_gaussian(:))/middle_contour_volume;
    mid_on_upper_raw_90_gaussian              = upper_contour_raw_90_gaussian.*middle_contour;
    lower_on_mid_raw_90_gaussian              = middle_contour.*lower_contour_raw_90_gaussian;
    upper_subset_mid_raw_90_gaussian          = upper_contour_raw_90_gaussian - mid_on_upper_raw_90_gaussian;
    mid_subset_lower_raw_90_gaussian          = middle_contour - lower_on_mid_raw_90_gaussian;    
    
    lower_contour_raw_95_gaussian             = observed_mean >= thr - supGa_raw_95_gaussian*tau*observed_std;
    upper_contour_raw_95_gaussian             = observed_mean >= thr + supGa_raw_95_gaussian*tau*observed_std;
    lower_contour_raw_95_gaussian_volume_prct = sum(lower_contour_raw_95_gaussian(:))/middle_contour_volume;
    upper_contour_raw_95_gaussian_volume_prct = sum(upper_contour_raw_95_gaussian(:))/middle_contour_volume;
    mid_on_upper_raw_95_gaussian              = upper_contour_raw_95_gaussian.*middle_contour;
    lower_on_mid_raw_95_gaussian              = middle_contour.*lower_contour_raw_95_gaussian;
    upper_subset_mid_raw_95_gaussian          = upper_contour_raw_95_gaussian - mid_on_upper_raw_95_gaussian;
    mid_subset_lower_raw_95_gaussian          = middle_contour - lower_on_mid_raw_95_gaussian;

    % Estimated Boundary
    supGa_raw_80_ero_dil_gaussian                     = prctile(supG_raw_ero_dil_gaussian, 80);
    supGa_raw_90_ero_dil_gaussian                     = prctile(supG_raw_ero_dil_gaussian, 90);
    supGa_raw_95_ero_dil_gaussian                     = prctile(supG_raw_ero_dil_gaussian, 95);
       
    lower_contour_raw_80_ero_dil_gaussian             = observed_mean >= thr - supGa_raw_80_ero_dil_gaussian*tau*observed_std;
    upper_contour_raw_80_ero_dil_gaussian             = observed_mean >= thr + supGa_raw_80_ero_dil_gaussian*tau*observed_std;
    lower_contour_raw_80_ero_dil_gaussian_volume_prct = sum(lower_contour_raw_80_ero_dil_gaussian(:))/middle_contour_volume;
    upper_contour_raw_80_ero_dil_gaussian_volume_prct = sum(upper_contour_raw_80_ero_dil_gaussian(:))/middle_contour_volume;
    mid_on_upper_raw_80_ero_dil_gaussian              = upper_contour_raw_80.*middle_contour;
    lower_on_mid_raw_80_ero_dil_gaussian              = middle_contour.*lower_contour_raw_80_ero_dil_gaussian;
    upper_subset_mid_raw_80_ero_dil_gaussian          = upper_contour_raw_80_ero_dil_gaussian - mid_on_upper_raw_80_ero_dil_gaussian;
    mid_subset_lower_raw_80_ero_dil_gaussian          = middle_contour - lower_on_mid_raw_80_ero_dil_gaussian;
    
    lower_contour_raw_90_ero_dil_gaussian             = observed_mean >= thr - supGa_raw_90_ero_dil_gaussian*tau*observed_std;
    upper_contour_raw_90_ero_dil_gaussian             = observed_mean >= thr + supGa_raw_90_ero_dil_gaussian*tau*observed_std;
    lower_contour_raw_90_ero_dil_gaussian_volume_prct = sum(lower_contour_raw_90_ero_dil_gaussian(:))/middle_contour_volume;
    upper_contour_raw_90_ero_dil_gaussian_volume_prct = sum(upper_contour_raw_90_ero_dil_gaussian(:))/middle_contour_volume;
    mid_on_upper_raw_90_ero_dil_gaussian              = upper_contour_raw_90_ero_dil_gaussian.*middle_contour;
    lower_on_mid_raw_90_ero_dil_gaussian              = middle_contour.*lower_contour_raw_90_ero_dil_gaussian;
    upper_subset_mid_raw_90_ero_dil_gaussian          = upper_contour_raw_90_ero_dil_gaussian - mid_on_upper_raw_90_ero_dil_gaussian;
    mid_subset_lower_raw_90_ero_dil_gaussian          = middle_contour - lower_on_mid_raw_90_ero_dil_gaussian;    
    
    lower_contour_raw_95_ero_dil_gaussian             = observed_mean >= thr - supGa_raw_95_ero_dil_gaussian*tau*observed_std;
    upper_contour_raw_95_ero_dil_gaussian             = observed_mean >= thr + supGa_raw_95_ero_dil_gaussian*tau*observed_std;
    lower_contour_raw_95_ero_dil_gaussian_volume_prct = sum(lower_contour_raw_95_ero_dil_gaussian(:))/middle_contour_volume;
    upper_contour_raw_95_ero_dil_gaussian_volume_prct = sum(upper_contour_raw_95_ero_dil_gaussian(:))/middle_contour_volume;
    mid_on_upper_raw_95_ero_dil_gaussian              = upper_contour_raw_95_ero_dil_gaussian.*middle_contour;
    lower_on_mid_raw_95_ero_dil_gaussian              = middle_contour.*lower_contour_raw_95_ero_dil_gaussian;
    upper_subset_mid_raw_95_ero_dil_gaussian          = upper_contour_raw_95_ero_dil_gaussian - mid_on_upper_raw_95_ero_dil_gaussian;
    mid_subset_lower_raw_95_ero_dil_gaussian          = middle_contour - lower_on_mid_raw_95_ero_dil_gaussian;

    %% Rademacher random variable results for the true and estimated boundary
    % True Boundary
    supGa_raw_80_signflips                     = prctile(supG_raw_signflips, 80);
    supGa_raw_90_signflips                     = prctile(supG_raw_signflips, 90);
    supGa_raw_95_signflips                     = prctile(supG_raw_signflips, 95);
       
    lower_contour_raw_80_signflips             = observed_mean >= thr - supGa_raw_80_signflips*tau*observed_std;
    upper_contour_raw_80_signflips             = observed_mean >= thr + supGa_raw_80_signflips*tau*observed_std;
    lower_contour_raw_80_signflips_volume_prct = sum(lower_contour_raw_80_signflips(:))/middle_contour_volume;
    upper_contour_raw_80_signflips_volume_prct = sum(upper_contour_raw_80_signflips(:))/middle_contour_volume;
    mid_on_upper_raw_80_signflips              = upper_contour_raw_80.*middle_contour;
    lower_on_mid_raw_80_signflips              = middle_contour.*lower_contour_raw_80_signflips;
    upper_subset_mid_raw_80_signflips          = upper_contour_raw_80_signflips - mid_on_upper_raw_80_signflips;
    mid_subset_lower_raw_80_signflips          = middle_contour - lower_on_mid_raw_80_signflips;
    
    lower_contour_raw_90_signflips             = observed_mean >= thr - supGa_raw_90_signflips*tau*observed_std;
    upper_contour_raw_90_signflips             = observed_mean >= thr + supGa_raw_90_signflips*tau*observed_std;
    lower_contour_raw_90_signflips_volume_prct = sum(lower_contour_raw_90_signflips(:))/middle_contour_volume;
    upper_contour_raw_90_signflips_volume_prct = sum(upper_contour_raw_90_signflips(:))/middle_contour_volume;
    mid_on_upper_raw_90_signflips              = upper_contour_raw_90_signflips.*middle_contour;
    lower_on_mid_raw_90_signflips              = middle_contour.*lower_contour_raw_90_signflips;
    upper_subset_mid_raw_90_signflips          = upper_contour_raw_90_signflips - mid_on_upper_raw_90_signflips;
    mid_subset_lower_raw_90_signflips          = middle_contour - lower_on_mid_raw_90_signflips;    
    
    lower_contour_raw_95_signflips             = observed_mean >= thr - supGa_raw_95_signflips*tau*observed_std;
    upper_contour_raw_95_signflips             = observed_mean >= thr + supGa_raw_95_signflips*tau*observed_std;
    lower_contour_raw_95_signflips_volume_prct = sum(lower_contour_raw_95_signflips(:))/middle_contour_volume;
    upper_contour_raw_95_signflips_volume_prct = sum(upper_contour_raw_95_signflips(:))/middle_contour_volume;
    mid_on_upper_raw_95_signflips              = upper_contour_raw_95_signflips.*middle_contour;
    lower_on_mid_raw_95_signflips              = middle_contour.*lower_contour_raw_95_signflips;
    upper_subset_mid_raw_95_signflips          = upper_contour_raw_95_signflips - mid_on_upper_raw_95_signflips;
    mid_subset_lower_raw_95_signflips          = middle_contour - lower_on_mid_raw_95_signflips;

    % Estimated Boundary
    supGa_raw_80_ero_dil_signflips                     = prctile(supG_raw_ero_dil_signflips, 80);
    supGa_raw_90_ero_dil_signflips                     = prctile(supG_raw_ero_dil_signflips, 90);
    supGa_raw_95_ero_dil_signflips                     = prctile(supG_raw_ero_dil_signflips, 95);
       
    lower_contour_raw_80_ero_dil_signflips             = observed_mean >= thr - supGa_raw_80_ero_dil_signflips*tau*observed_std;
    upper_contour_raw_80_ero_dil_signflips             = observed_mean >= thr + supGa_raw_80_ero_dil_signflips*tau*observed_std;
    lower_contour_raw_80_ero_dil_signflips_volume_prct = sum(lower_contour_raw_80_ero_dil_signflips(:))/middle_contour_volume;
    upper_contour_raw_80_ero_dil_signflips_volume_prct = sum(upper_contour_raw_80_ero_dil_signflips(:))/middle_contour_volume;
    mid_on_upper_raw_80_ero_dil_signflips              = upper_contour_raw_80.*middle_contour;
    lower_on_mid_raw_80_ero_dil_signflips              = middle_contour.*lower_contour_raw_80_ero_dil_signflips;
    upper_subset_mid_raw_80_ero_dil_signflips          = upper_contour_raw_80_ero_dil_signflips - mid_on_upper_raw_80_ero_dil_signflips;
    mid_subset_lower_raw_80_ero_dil_signflips          = middle_contour - lower_on_mid_raw_80_ero_dil_signflips;
    
    lower_contour_raw_90_ero_dil_signflips             = observed_mean >= thr - supGa_raw_90_ero_dil_signflips*tau*observed_std;
    upper_contour_raw_90_ero_dil_signflips             = observed_mean >= thr + supGa_raw_90_ero_dil_signflips*tau*observed_std;
    lower_contour_raw_90_ero_dil_signflips_volume_prct = sum(lower_contour_raw_90_ero_dil_signflips(:))/middle_contour_volume;
    upper_contour_raw_90_ero_dil_signflips_volume_prct = sum(upper_contour_raw_90_ero_dil_signflips(:))/middle_contour_volume;
    mid_on_upper_raw_90_ero_dil_signflips              = upper_contour_raw_90_ero_dil_signflips.*middle_contour;
    lower_on_mid_raw_90_ero_dil_signflips              = middle_contour.*lower_contour_raw_90_ero_dil_signflips;
    upper_subset_mid_raw_90_ero_dil_signflips          = upper_contour_raw_90_ero_dil_signflips - mid_on_upper_raw_90_ero_dil_signflips;
    mid_subset_lower_raw_90_ero_dil_signflips          = middle_contour - lower_on_mid_raw_90_ero_dil_signflips;    
    
    lower_contour_raw_95_ero_dil_signflips             = observed_mean >= thr - supGa_raw_95_ero_dil_signflips*tau*observed_std;
    upper_contour_raw_95_ero_dil_signflips             = observed_mean >= thr + supGa_raw_95_ero_dil_signflips*tau*observed_std;
    lower_contour_raw_95_ero_dil_signflips_volume_prct = sum(lower_contour_raw_95_ero_dil_signflips(:))/middle_contour_volume;
    upper_contour_raw_95_ero_dil_signflips_volume_prct = sum(upper_contour_raw_95_ero_dil_signflips(:))/middle_contour_volume;
    mid_on_upper_raw_95_ero_dil_signflips              = upper_contour_raw_95_ero_dil_signflips.*middle_contour;
    lower_on_mid_raw_95_ero_dil_signflips              = middle_contour.*lower_contour_raw_95_ero_dil_signflips;
    upper_subset_mid_raw_95_ero_dil_signflips          = upper_contour_raw_95_ero_dil_signflips - mid_on_upper_raw_95_ero_dil_signflips;
    mid_subset_lower_raw_95_ero_dil_signflips          = middle_contour - lower_on_mid_raw_95_ero_dil_signflips;
    
    %
    % Storing all variables of interest
    %
    supG_raw_gaussian_store(:,t)                  = supG_raw_gaussian;
    supG_raw_ero_dil_gaussian_store(:,t)          = supG_raw_ero_dil_gaussian;
    supG_raw_signflips_store(:,t)                 = supG_raw_signflips;
    supG_raw_ero_dil_signflips_store(:,t)         = supG_raw_ero_dil_signflips;
      
    threshold_raw_80_gaussian_store(t)                              = supGa_raw_80_gaussian;
    lower_contour_raw_80_gaussian_store(t,:,:)                      = lower_contour_raw_80_gaussian;
    upper_contour_raw_80_gaussian_store(t,:,:)                      = upper_contour_raw_80_gaussian;
    upper_subset_mid_raw_80_gaussian_store(t,:,:)                   = upper_subset_mid_raw_80_gaussian;
    mid_subset_lower_raw_80_gaussian_store(t,:,:)                   = mid_subset_lower_raw_80_gaussian;
    lower_contour_raw_80_gaussian_volume_prct_store(t)              = lower_contour_raw_80_gaussian_volume_prct;
    upper_contour_raw_80_gaussian_volume_prct_store(t)              = upper_contour_raw_80_gaussian_volume_prct;
    threshold_raw_80_ero_dil_gaussian_store(t)                      = supGa_raw_80_ero_dil_gaussian;
    lower_contour_raw_80_ero_dil_gaussian_store(t,:,:)              = lower_contour_raw_80_ero_dil_gaussian;
    upper_contour_raw_80_ero_dil_gaussian_store(t,:,:)              = upper_contour_raw_80_ero_dil_gaussian;
    upper_subset_mid_raw_80_ero_dil_gaussian_store(t,:,:)           = upper_subset_mid_raw_80_ero_dil_gaussian;
    mid_subset_lower_raw_80_ero_dil_gaussian_store(t,:,:)           = mid_subset_lower_raw_80_ero_dil_gaussian;
    lower_contour_raw_80_ero_dil_gaussian_volume_prct_store(t)      = lower_contour_raw_80_ero_dil_gaussian_volume_prct;
    upper_contour_raw_80_ero_dil_gaussian_volume_prct_store(t)      = upper_contour_raw_80_ero_dil_gaussian_volume_prct;
    threshold_raw_80_signflips_store(t)                              = supGa_raw_80_signflips;
    lower_contour_raw_80_signflips_store(t,:,:)                      = lower_contour_raw_80_signflips;
    upper_contour_raw_80_signflips_store(t,:,:)                      = upper_contour_raw_80_signflips;
    upper_subset_mid_raw_80_signflips_store(t,:,:)                   = upper_subset_mid_raw_80_signflips;
    mid_subset_lower_raw_80_signflips_store(t,:,:)                   = mid_subset_lower_raw_80_signflips;
    lower_contour_raw_80_signflips_volume_prct_store(t)              = lower_contour_raw_80_signflips_volume_prct;
    upper_contour_raw_80_signflips_volume_prct_store(t)              = upper_contour_raw_80_signflips_volume_prct;
    threshold_raw_80_ero_dil_signflips_store(t)                      = supGa_raw_80_ero_dil_signflips;
    lower_contour_raw_80_ero_dil_signflips_store(t,:,:)              = lower_contour_raw_80_ero_dil_signflips;
    upper_contour_raw_80_ero_dil_signflips_store(t,:,:)              = upper_contour_raw_80_ero_dil_signflips;
    upper_subset_mid_raw_80_ero_dil_signflips_store(t,:,:)           = upper_subset_mid_raw_80_ero_dil_signflips;
    mid_subset_lower_raw_80_ero_dil_signflips_store(t,:,:)           = mid_subset_lower_raw_80_ero_dil_signflips;
    lower_contour_raw_80_ero_dil_signflips_volume_prct_store(t)      = lower_contour_raw_80_ero_dil_signflips_volume_prct;
    upper_contour_raw_80_ero_dil_signflips_volume_prct_store(t)      = upper_contour_raw_80_ero_dil_signflips_volume_prct;

    threshold_raw_90_gaussian_store(t)                              = supGa_raw_90_gaussian;
    lower_contour_raw_90_gaussian_store(t,:,:)                      = lower_contour_raw_90_gaussian;
    upper_contour_raw_90_gaussian_store(t,:,:)                      = upper_contour_raw_90_gaussian;
    upper_subset_mid_raw_90_gaussian_store(t,:,:)                   = upper_subset_mid_raw_90_gaussian;
    mid_subset_lower_raw_90_gaussian_store(t,:,:)                   = mid_subset_lower_raw_90_gaussian;
    lower_contour_raw_90_gaussian_volume_prct_store(t)              = lower_contour_raw_90_gaussian_volume_prct;
    upper_contour_raw_90_gaussian_volume_prct_store(t)              = upper_contour_raw_90_gaussian_volume_prct;
    threshold_raw_90_ero_dil_gaussian_store(t)                      = supGa_raw_90_ero_dil_gaussian;
    lower_contour_raw_90_ero_dil_gaussian_store(t,:,:)              = lower_contour_raw_90_ero_dil_gaussian;
    upper_contour_raw_90_ero_dil_gaussian_store(t,:,:)              = upper_contour_raw_90_ero_dil_gaussian;
    upper_subset_mid_raw_90_ero_dil_gaussian_store(t,:,:)           = upper_subset_mid_raw_90_ero_dil_gaussian;
    mid_subset_lower_raw_90_ero_dil_gaussian_store(t,:,:)           = mid_subset_lower_raw_90_ero_dil_gaussian;
    lower_contour_raw_90_ero_dil_gaussian_volume_prct_store(t)      = lower_contour_raw_90_ero_dil_gaussian_volume_prct;
    upper_contour_raw_90_ero_dil_gaussian_volume_prct_store(t)      = upper_contour_raw_90_ero_dil_gaussian_volume_prct;
    threshold_raw_90_signflips_store(t)                              = supGa_raw_90_signflips;
    lower_contour_raw_90_signflips_store(t,:,:)                      = lower_contour_raw_90_signflips;
    upper_contour_raw_90_signflips_store(t,:,:)                      = upper_contour_raw_90_signflips;
    upper_subset_mid_raw_90_signflips_store(t,:,:)                   = upper_subset_mid_raw_90_signflips;
    mid_subset_lower_raw_90_signflips_store(t,:,:)                   = mid_subset_lower_raw_90_signflips;
    lower_contour_raw_90_signflips_volume_prct_store(t)              = lower_contour_raw_90_signflips_volume_prct;
    upper_contour_raw_90_signflips_volume_prct_store(t)              = upper_contour_raw_90_signflips_volume_prct;
    threshold_raw_90_ero_dil_signflips_store(t)                      = supGa_raw_90_ero_dil_signflips;
    lower_contour_raw_90_ero_dil_signflips_store(t,:,:)              = lower_contour_raw_90_ero_dil_signflips;
    upper_contour_raw_90_ero_dil_signflips_store(t,:,:)              = upper_contour_raw_90_ero_dil_signflips;
    upper_subset_mid_raw_90_ero_dil_signflips_store(t,:,:)           = upper_subset_mid_raw_90_ero_dil_signflips;
    mid_subset_lower_raw_90_ero_dil_signflips_store(t,:,:)           = mid_subset_lower_raw_90_ero_dil_signflips;
    lower_contour_raw_90_ero_dil_signflips_volume_prct_store(t)      = lower_contour_raw_90_ero_dil_signflips_volume_prct;
    upper_contour_raw_90_ero_dil_signflips_volume_prct_store(t)      = upper_contour_raw_90_ero_dil_signflips_volume_prct;
    
    threshold_raw_95_gaussian_store(t)                              = supGa_raw_95_gaussian;
    lower_contour_raw_95_gaussian_store(t,:,:)                      = lower_contour_raw_95_gaussian;
    upper_contour_raw_95_gaussian_store(t,:,:)                      = upper_contour_raw_95_gaussian;
    upper_subset_mid_raw_95_gaussian_store(t,:,:)                   = upper_subset_mid_raw_95_gaussian;
    mid_subset_lower_raw_95_gaussian_store(t,:,:)                   = mid_subset_lower_raw_95_gaussian;
    threshold_raw_95_ero_dil_gaussian_store(t)                      = supGa_raw_95_ero_dil_gaussian;
    lower_contour_raw_95_ero_dil_gaussian_store(t,:,:)              = lower_contour_raw_95_ero_dil_gaussian;
    upper_contour_raw_95_ero_dil_gaussian_store(t,:,:)              = upper_contour_raw_95_ero_dil_gaussian;
    upper_subset_mid_raw_95_ero_dil_gaussian_store(t,:,:)           = upper_subset_mid_raw_95_ero_dil_gaussian;
    mid_subset_lower_raw_95_ero_dil_gaussian_store(t,:,:)           = mid_subset_lower_raw_95_ero_dil_gaussian;
    threshold_raw_95_signflips_store(t)                              = supGa_raw_95_signflips;
    lower_contour_raw_95_signflips_store(t,:,:)                      = lower_contour_raw_95_signflips;
    upper_contour_raw_95_signflips_store(t,:,:)                      = upper_contour_raw_95_signflips;
    upper_subset_mid_raw_95_signflips_store(t,:,:)                   = upper_subset_mid_raw_95_signflips;
    mid_subset_lower_raw_95_signflips_store(t,:,:)                   = mid_subset_lower_raw_95_signflips;
    threshold_raw_95_ero_dil_signflips_store(t)                      = supGa_raw_95_ero_dil_signflips;
    lower_contour_raw_95_ero_dil_signflips_store(t,:,:)              = lower_contour_raw_95_ero_dil_signflips;
    upper_contour_raw_95_ero_dil_signflips_store(t,:,:)              = upper_contour_raw_95_ero_dil_signflips;
    upper_subset_mid_raw_95_ero_dil_signflips_store(t,:,:)           = upper_subset_mid_raw_95_ero_dil_signflips;
    mid_subset_lower_raw_95_ero_dil_signflips_store(t,:,:)           = mid_subset_lower_raw_95_ero_dil_signflips;

    threshold_raw_95_gaussian_store(t)                              = supGa_raw_95_gaussian;
    lower_contour_raw_95_gaussian_store(t,:,:)                      = lower_contour_raw_95_gaussian;
    upper_contour_raw_95_gaussian_store(t,:,:)                      = upper_contour_raw_95_gaussian;
    upper_subset_mid_raw_95_gaussian_store(t,:,:)                   = upper_subset_mid_raw_95_gaussian;
    mid_subset_lower_raw_95_gaussian_store(t,:,:)                   = mid_subset_lower_raw_95_gaussian;
    lower_contour_raw_95_gaussian_volume_prct_store(t)              = lower_contour_raw_95_gaussian_volume_prct;
    upper_contour_raw_95_gaussian_volume_prct_store(t)              = upper_contour_raw_95_gaussian_volume_prct;
    threshold_raw_95_ero_dil_gaussian_store(t)                      = supGa_raw_95_ero_dil_gaussian;
    lower_contour_raw_95_ero_dil_gaussian_store(t,:,:)              = lower_contour_raw_95_ero_dil_gaussian;
    upper_contour_raw_95_ero_dil_gaussian_store(t,:,:)              = upper_contour_raw_95_ero_dil_gaussian;
    upper_subset_mid_raw_95_ero_dil_gaussian_store(t,:,:)           = upper_subset_mid_raw_95_ero_dil_gaussian;
    mid_subset_lower_raw_95_ero_dil_gaussian_store(t,:,:)           = mid_subset_lower_raw_95_ero_dil_gaussian;
    lower_contour_raw_95_ero_dil_gaussian_volume_prct_store(t)      = lower_contour_raw_95_ero_dil_gaussian_volume_prct;
    upper_contour_raw_95_ero_dil_gaussian_volume_prct_store(t)      = upper_contour_raw_95_ero_dil_gaussian_volume_prct;
    threshold_raw_95_signflips_store(t)                              = supGa_raw_95_signflips;
    lower_contour_raw_95_signflips_store(t,:,:)                      = lower_contour_raw_95_signflips;
    upper_contour_raw_95_signflips_store(t,:,:)                      = upper_contour_raw_95_signflips;
    upper_subset_mid_raw_95_signflips_store(t,:,:)                   = upper_subset_mid_raw_95_signflips;
    mid_subset_lower_raw_95_signflips_store(t,:,:)                   = mid_subset_lower_raw_95_signflips;
    lower_contour_raw_95_signflips_volume_prct_store(t)              = lower_contour_raw_95_signflips_volume_prct;
    upper_contour_raw_95_signflips_volume_prct_store(t)              = upper_contour_raw_95_signflips_volume_prct;
    threshold_raw_95_ero_dil_signflips_store(t)                      = supGa_raw_95_ero_dil_signflips;
    lower_contour_raw_95_ero_dil_signflips_store(t,:,:)              = lower_contour_raw_95_ero_dil_signflips;
    upper_contour_raw_95_ero_dil_signflips_store(t,:,:)              = upper_contour_raw_95_ero_dil_signflips;
    upper_subset_mid_raw_95_ero_dil_signflips_store(t,:,:)           = upper_subset_mid_raw_95_ero_dil_signflips;
    mid_subset_lower_raw_95_ero_dil_signflips_store(t,:,:)           = mid_subset_lower_raw_95_ero_dil_signflips;
    lower_contour_raw_95_ero_dil_signflips_volume_prct_store(t)      = lower_contour_raw_95_ero_dil_signflips_volume_prct;
    upper_contour_raw_95_ero_dil_signflips_volume_prct_store(t)      = upper_contour_raw_95_ero_dil_signflips_volume_prct;
    
    if sum(upper_subset_mid_raw_80_gaussian(:))+sum(mid_subset_lower_raw_80_gaussian(:))==0
      subset_success_vector_raw_80_gaussian(t) = 1; 
      fprintf('raw nominal 80 gaussian true boundary success! \n');
    else 
      subset_success_vector_raw_80_gaussian(t) = 0; 
      fprintf('raw nominal 80 gaussian true boundary failure! \n');
    end 
    
    if sum(upper_subset_mid_raw_80_ero_dil_gaussian(:))+sum(mid_subset_lower_raw_80_ero_dil_gaussian(:))==0
      subset_success_vector_raw_80_ero_dil_gaussian(t) = 1; 
      fprintf('raw nominal 80 gaussian estimated boundary success! \n');
    else 
      subset_success_vector_raw_80_ero_dil_gaussian(t) = 0; 
      fprintf('raw nominal 80 gaussian estimated boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_80_signflips(:))+sum(mid_subset_lower_raw_80_signflips(:))==0
      subset_success_vector_raw_80_signflips(t) = 1; 
      fprintf('raw nominal 80 signflips estimated boundary success! \n');
    else 
      subset_success_vector_raw_80_signflips(t) = 0; 
      fprintf('raw nominal 80 signflips estimated boundary failure! \n');
    end 
    
    if sum(upper_subset_mid_raw_80_ero_dil_signflips(:))+sum(mid_subset_lower_raw_80_ero_dil_signflips(:))==0
      subset_success_vector_raw_80_ero_dil_signflips(t) = 1; 
      fprintf('raw nominal 80 signflips estimated boundary success! \n');
    else 
      subset_success_vector_raw_80_ero_dil_signflips(t) = 0; 
      fprintf('raw nominal 80 signflips estimated boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_90_gaussian(:))+sum(mid_subset_lower_raw_90_gaussian(:))==0
      subset_success_vector_raw_90_gaussian(t) = 1; 
      fprintf('raw nominal 90 gaussian true boundary success! \n');
    else 
      subset_success_vector_raw_90_gaussian(t) = 0; 
      fprintf('raw nominal 90 gaussian true boundary failure! \n');
    end 
    
    if sum(upper_subset_mid_raw_90_ero_dil_gaussian(:))+sum(mid_subset_lower_raw_90_ero_dil_gaussian(:))==0
      subset_success_vector_raw_90_ero_dil_gaussian(t) = 1; 
      fprintf('raw nominal 90 gaussian estimated boundary success! \n');
    else 
      subset_success_vector_raw_90_ero_dil_gaussian(t) = 0; 
      fprintf('raw nominal 90 gaussian estimated boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_90_signflips(:))+sum(mid_subset_lower_raw_90_signflips(:))==0
      subset_success_vector_raw_90_signflips(t) = 1; 
      fprintf('raw nominal 90 signflips estimated boundary success! \n');
    else 
      subset_success_vector_raw_90_signflips(t) = 0; 
      fprintf('raw nominal 90 signflips estimated boundary failure! \n');
    end 
    
    if sum(upper_subset_mid_raw_90_ero_dil_signflips(:))+sum(mid_subset_lower_raw_90_ero_dil_signflips(:))==0
      subset_success_vector_raw_90_ero_dil_signflips(t) = 1; 
      fprintf('raw nominal 90 signflips estimated boundary success! \n');
    else 
      subset_success_vector_raw_90_ero_dil_signflips(t) = 0; 
      fprintf('raw nominal 90 signflips estimated boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_95_gaussian(:))+sum(mid_subset_lower_raw_95_gaussian(:))==0
      subset_success_vector_raw_95_gaussian(t) = 1; 
      fprintf('raw nominal 95 gaussian true boundary success! \n');
    else 
      subset_success_vector_raw_95_gaussian(t) = 0; 
      fprintf('raw nominal 95 gaussian true boundary failure! \n');
    end 
    
    if sum(upper_subset_mid_raw_95_ero_dil_gaussian(:))+sum(mid_subset_lower_raw_95_ero_dil_gaussian(:))==0
      subset_success_vector_raw_95_ero_dil_gaussian(t) = 1; 
      fprintf('raw nominal 95 gaussian estimated boundary success! \n');
    else 
      subset_success_vector_raw_95_ero_dil_gaussian(t) = 0; 
      fprintf('raw nominal 95 gaussian estimated boundary failure! \n');
    end 

    if sum(upper_subset_mid_raw_95_signflips(:))+sum(mid_subset_lower_raw_95_signflips(:))==0
      subset_success_vector_raw_95_signflips(t) = 1; 
      fprintf('raw nominal 95 signflips estimated boundary success! \n');
    else 
      subset_success_vector_raw_95_signflips(t) = 0; 
      fprintf('raw nominal 95 signflips estimated boundary failure! \n');
    end 
    
    if sum(upper_subset_mid_raw_95_ero_dil_signflips(:))+sum(mid_subset_lower_raw_95_ero_dil_signflips(:))==0
      subset_success_vector_raw_95_ero_dil_signflips(t) = 1; 
      fprintf('raw nominal 95 signflips estimated boundary success! \n');
    else 
      subset_success_vector_raw_95_ero_dil_signflips(t) = 0; 
      fprintf('raw nominal 95 signflips estimated boundary failure! \n');
    end     
end

percentage_success_vector_raw_80_gaussian                    = mean(subset_success_vector_raw_80_gaussian, 1);
percentage_success_vector_raw_80_ero_dil_gaussian            = mean(subset_success_vector_raw_80_ero_dil_gaussian, 1);
percentage_success_vector_raw_80_signflips                   = mean(subset_success_vector_raw_80_signflips, 1);
percentage_success_vector_raw_80_ero_dil_signflips           = mean(subset_success_vector_raw_80_ero_dil_signflips, 1);

percentage_success_vector_raw_90_gaussian                    = mean(subset_success_vector_raw_90_gaussian, 1);
percentage_success_vector_raw_90_ero_dil_gaussian            = mean(subset_success_vector_raw_90_ero_dil_gaussian, 1);
percentage_success_vector_raw_90_signflips                   = mean(subset_success_vector_raw_90_signflips, 1);
percentage_success_vector_raw_90_ero_dil_signflips           = mean(subset_success_vector_raw_90_ero_dil_signflips, 1);

percentage_success_vector_raw_95_gaussian                    = mean(subset_success_vector_raw_95_gaussian, 1);
percentage_success_vector_raw_95_ero_dil_gaussian            = mean(subset_success_vector_raw_95_ero_dil_gaussian, 1);
percentage_success_vector_raw_95_signflips                   = mean(subset_success_vector_raw_95_signflips, 1);
percentage_success_vector_raw_95_ero_dil_signflips           = mean(subset_success_vector_raw_95_ero_dil_signflips, 1);


eval(['save ' SvNm ' nSubj nRlz dim smo mag rimFWHM thr nBoot '... 
      'threshold_raw_80_gaussian_store threshold_raw_90_gaussian_store threshold_raw_95_gaussian_store threshold_raw_80_signflips_store threshold_raw_90_signflips_store threshold_raw_95_signflips_store threshold_raw_80_ero_dil_gaussian_store threshold_raw_90_ero_dil_gaussian_store threshold_raw_95_ero_dil_gaussian_store threshold_raw_80_ero_dil_signflips_store threshold_raw_90_ero_dil_signflips_store threshold_raw_95_ero_dil_signflips_store '...
      'lower_contour_raw_80_gaussian_store lower_contour_raw_90_gaussian_store lower_contour_raw_95_gaussian_store lower_contour_raw_80_signflips_store lower_contour_raw_90_signflips_store lower_contour_raw_95_signflips_store lower_contour_raw_80_ero_dil_gaussian_store lower_contour_raw_90_ero_dil_gaussian_store lower_contour_raw_95_ero_dil_gaussian_store lower_contour_raw_80_ero_dil_signflips_store lower_contour_raw_90_ero_dil_signflips_store lower_contour_raw_95_ero_dil_signflips_store '...
      'upper_contour_raw_80_gaussian_store upper_contour_raw_90_gaussian_store upper_contour_raw_95_gaussian_store upper_contour_raw_80_signflips_store upper_contour_raw_90_signflips_store upper_contour_raw_95_signflips_store upper_contour_raw_80_ero_dil_gaussian_store upper_contour_raw_90_ero_dil_gaussian_store upper_contour_raw_95_ero_dil_gaussian_store upper_contour_raw_80_ero_dil_signflips_store upper_contour_raw_90_ero_dil_signflips_store upper_contour_raw_95_ero_dil_signflips_store '...
      'upper_subset_mid_raw_80_gaussian_store upper_subset_mid_raw_90_gaussian_store upper_subset_mid_raw_95_gaussian_store upper_subset_mid_raw_80_signflips_store upper_subset_mid_raw_90_signflips_store upper_subset_mid_raw_95_signflips_store upper_subset_mid_raw_80_ero_dil_gaussian_store upper_subset_mid_raw_90_ero_dil_gaussian_store upper_subset_mid_raw_95_ero_dil_gaussian_store upper_subset_mid_raw_80_ero_dil_signflips_store upper_subset_mid_raw_90_ero_dil_signflips_store upper_subset_mid_raw_95_ero_dil_signflips_store '...
      'mid_subset_lower_raw_80_gaussian_store mid_subset_lower_raw_90_gaussian_store mid_subset_lower_raw_95_gaussian_store mid_subset_lower_raw_80_signflips_store mid_subset_lower_raw_90_signflips_store mid_subset_lower_raw_95_signflips_store mid_subset_lower_raw_80_ero_dil_gaussian_store mid_subset_lower_raw_90_ero_dil_gaussian_store mid_subset_lower_raw_95_ero_dil_gaussian_store mid_subset_lower_raw_80_ero_dil_signflips_store mid_subset_lower_raw_90_ero_dil_signflips_store mid_subset_lower_raw_95_ero_dil_signflips_store '...
      'subset_success_vector_raw_80_gaussian subset_success_vector_raw_90_gaussian subset_success_vector_raw_95_gaussian subset_success_vector_raw_80_signflips subset_success_vector_raw_90_signflips subset_success_vector_raw_95_signflips subset_success_vector_raw_80_ero_dil_gaussian subset_success_vector_raw_90_ero_dil_gaussian subset_success_vector_raw_95_ero_dil_gaussian subset_success_vector_raw_80_ero_dil_signflips subset_success_vector_raw_90_ero_dil_signflips subset_success_vector_raw_95_ero_dil_signflips '...
      'percentage_success_vector_raw_80_gaussian percentage_success_vector_raw_90_gaussian percentage_success_vector_raw_95_gaussian percentage_success_vector_raw_80_signflips percentage_success_vector_raw_90_signflips percentage_success_vector_raw_95_signflips percentage_success_vector_raw_80_ero_dil_gaussian percentage_success_vector_raw_90_ero_dil_gaussian percentage_success_vector_raw_95_ero_dil_gaussian percentage_success_vector_raw_80_ero_dil_signflips percentage_success_vector_raw_90_ero_dil_signflips percentage_success_vector_raw_95_ero_dil_signflips '...
      'supG_raw_gaussian_store supG_raw_signflips_store supG_raw_ero_dil_gaussian_store supG_raw_ero_dil_signflips '...
      'raw_field_boundary_store '...
      'middle_contour_volume observed_AC_volume '...
      'lower_contour_raw_80_gaussian_volume_prct_store lower_contour_raw_90_gaussian_volume_prct_store lower_contour_raw_95_gaussian_volume_prct_store lower_contour_raw_80_ero_dil_gaussian_volume_prct_store lower_contour_raw_90_ero_dil_gaussian_volume_prct_store lower_contour_raw_95_ero_dil_gaussian_volume_prct_store lower_contour_raw_80_signflips_volume_prct_store lower_contour_raw_90_signflips_volume_prct_store lower_contour_raw_95_gaussian_volume_prct_store lower_contour_raw_80_ero_dil_signflips_volume_prct_store lower_contour_raw_90_ero_dil_signflips_volume_prct_store lower_contour_raw_95_ero_dil_signflips_volume_prct_store '...
      'upper_contour_raw_80_gaussian_volume_prct_store upper_contour_raw_90_gaussian_volume_prct_store upper_contour_raw_95_gaussian_volume_prct_store upper_contour_raw_80_ero_dil_gaussian_volume_prct_store upper_contour_raw_90_ero_dil_gaussian_volume_prct_store upper_contour_raw_95_ero_dil_gaussian_volume_prct_store upper_contour_raw_80_signflips_volume_prct_store upper_contour_raw_90_signflips_volume_prct_store upper_contour_raw_95_gaussian_volume_prct_store upper_contour_raw_80_ero_dil_signflips_volume_prct_store upper_contour_raw_90_ero_dil_signflips_volume_prct_store upper_contour_raw_95_ero_dil_signflips_volume_prct_store'])
