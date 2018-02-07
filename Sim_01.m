function Sim_01(nSubj,SvNm,nRlz)
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
nBoot   = 1000;
dim     = [100 100]; 
mag     = 4;
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
raw_noise    = zeros([wdim nSubj]);

% This stores the vector SupG for each run
% This vector stores the result for each realisation on whether AC^+ < AC < AC^ for each level of smoothing (1 if true, 0 if false) 
subset_success_vector_raw   = zeros(nRlz, 1); 
subset_success_vector_cohen   = zeros(nRlz, 1);

%- This vector stores the threshold value 'c' for each run
threshold_raw_store          = zeros(nRlz, 1);
threshold_cohen_store          = zeros(nRlz, 1);

% This stores the vector SupG for each run
supG_raw_store        = zeros(nBoot, nRlz);
supG_cohen_store        = zeros(nBoot, nRlz);

%-These matrices store all the sets of interest during the bootstrap
% method for all levels of smoothing
lower_contour_raw_store       = zeros([nRlz dim]);
upper_contour_raw_store       = zeros([nRlz dim]);
upper_subset_mid_raw_store = zeros([nRlz dim]);
mid_subset_lower_raw_store = zeros([nRlz dim]);

lower_contour_cohen_store       = zeros([nRlz dim]);
upper_contour_cohen_store       = zeros([nRlz dim]);
upper_subset_mid_cohen_store = zeros([nRlz dim]);
mid_subset_lower_cohen_store = zeros([nRlz dim]);

supG_raw         = zeros(nBoot,1);
supG_cohen         = zeros(nBoot,1);

raw_field_boundary_store        = zeros(dim(2), nBoot);
snr_field_boundary_store        = zeros(dim(2), nBoot);

% Creating linearly increasing signal across columns
Sig = repmat(linspace(0, mag), dim(2), 1);

% Uncomment to look at the Signal
% imagesc(Sig); axis image; colorbar
AC = Sig >= thr;

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
         mystd(reshape(observed_data,[prod(dim) nSubj]),stdblk),...
           dim); 
      
      observed_d = observed_mean./observed_std;
 
      % Residuals
      resid = bsxfun(@minus,observed_data,observed_mean);
      resid = spdiags(1./reshape(observed_std, [prod(dim) 1]), 0,prod(dim),prod(dim))*reshape(resid,[prod(dim) nSubj]);
      
      % SNR Residuals for implementation of Fabians method
      snr_resid = resid - 0.5*bsxfun(@times, resid.^2 - 1,reshape(observed_d, [prod(dim) 1]));
      snr_resid = spdiags(1./sqrt(1 + 0.5*observed_d.^2), 0, prod(dim), prod(dim))*snr_resid; 
      
      % Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          signflips = randi(2, [nSubj,1])*2-3;
          resid_bootstrap      = resid*spdiags(signflips, nSubj, nSubj);
          resid_bootstrap      = reshape(resid_bootstrap, [dim nSubj]);
          resid_field          = sum(resid_bootstrap, 3)/sqrt(nSubj); 
          
          snr_resid_bootstrap  = snr_resid*spdiags(signflips, nSubj, nSubj);
          snr_resid_bootstrap  = reshape(snr_resid_bootstrap, [dim nSubj]);
          snr_resid_field      = sum(snr_resid_bootstrap, 3)/sqrt(nSubj);

          supG_raw(k)        = max(resid_field(true_boundary));
          supG_cohen(k)        = min(snr_resid_field(true_boundary));
          
          if t==nRlz
              raw_field_boundary_store(:,k)         = resid_field(true_boundary);
              snr_field_boundary_store(:,k)         = snr_resid_field(true_boundary);
          end 
      end
      
    supGa_raw                     = prctile(supGa_raw, 95);
    supGa_cohen                   = prctile(supGa_snr, 95);
    
    middle_contour                = AC;
    
    lower_contour_raw             = observed_mean >= thr - supGa_raw*tau*observed_std;
    upper_contour_raw             = observed_mean >= thr + supGa_raw*tau*observed_std;
    mid_on_upper_raw              = upper_contour_raw.*middle_contour;
    lower_on_mid_raw              = middle_contour.*lower_contour_raw;
    upper_subset_mid_raw          = upper_contour_raw - mid_on_upper_raw;
    mid_subset_lower_raw          = middle_contour - lower_on_mid_raw;
    
    lower_contour_cohen           = observed_d >= thr - supGa_cohen*tau*sqrt(1 + 0.5*observed_d.^2);
    upper_contour_cohen           = observed_d >= thr + supGa_cohen*tau*sqrt(1 + 0.5*observed_d.^2);
    mid_on_upper_cohen            = upper_contour_cohen.*middle_contour;
    lower_on_mid_cohen            = middle_contour.*lower_contour_cohen;
    upper_subset_mid_cohen        = upper_contour_cohen - mid_on_upper_cohen;
    mid_subset_lower_cohen        = middle_contour - lower_on_mid_cohen;
    
    %
    % Storing all variables of interest
    %
    supG_raw_store(:,t)                 = supG_raw;
    threshold_raw_store(t)              = supGa_raw;
    lower_contour_raw_store(t,:,:)      = lower_contour_raw;
    upper_contour_raw_store(t,:,:)      = upper_contour_raw;
    upper_subset_mid_raw_store(t,:,:)   = upper_subset_mid_raw;
    mid_subset_lower_raw_store(t,:,:)   = mid_subset_lower_raw;
    
    supG_cohen_store(:,t)               = supG_cohen;
    threshold_cohen_store(t)            = supGa_cohen;
    lower_contour_cohen_store(t,:,:)    = lower_contour_cohen;
    upper_contour_cohen_store(t,:,:)    = upper_contour_cohen;
    upper_subset_mid_cohen_store(t,:,:) = upper_subset_mid_cohen;
    mid_subset_lower_cohen_store(t,:,:) = mid_subset_lower_cohen;
    
    if sum(upper_subset_mid_raw(:))+sum(mid_subset_lower_raw(:))==0
      subset_success_vector_raw(t) = 1; 
      fprintf('raw success! \n');
    else 
      subset_success_vector_raw(t) = 0; 
      fprintf('raw failure! \n');
    end 
    
    if sum(upper_subset_mid_cohen(:))+sum(mid_subset_lower_cohen(:))==0
      subset_success_vector_cohen(t) = 1; 
      fprintf('cohen success! \n');
    else 
      subset_success_vector_cohen(t) = 0; 
      fprintf('cohen failure! \n');
    end 
    
end

percentage_success_vector_raw         = mean(subset_success_vector_raw, 1);
percentage_success_vector_cohen       = mean(subset_success_vector_cohen, 1);

eval(['save ' SvNm ' nSubj nRlz dim smo mag rimFWHM thr nBoot threshold_raw_store_ threshold_cohen_store lower_contour_raw_store lower_contour_cohen_store upper_contour_raw_store upper_contour_cohen_store upper_subset_mid_raw_store upper_subset_mid_cohen_store mid_subset_lower_raw_store mid_subset_lower_cohen_store subset_success_vector_raw subset_success_vector_cohen percentage_success_vector_raw percentage_success_vector_cohen supG_raw_store supG_cohen_store raw_field_boundary_store snr_field_boundary_store'])