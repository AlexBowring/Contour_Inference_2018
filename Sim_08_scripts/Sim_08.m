function Sim_08(nSubj,SvNm,nRlz)
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
subset_success_vector_raw_80   = zeros(nRlz, 1); 
subset_success_vector_raw_90   = zeros(nRlz, 1);
subset_success_vector_raw_95   = zeros(nRlz, 1); 

subset_success_vector_raw_80_boundary   = zeros(nRlz, 1); 
subset_success_vector_raw_90_boundary   = zeros(nRlz, 1);
subset_success_vector_raw_95_boundary   = zeros(nRlz, 1); 

%- This vector stores the threshold value 'c' for each run
threshold_raw_80_store          = zeros(nRlz, 1);
threshold_raw_90_store          = zeros(nRlz, 1);
threshold_raw_95_store          = zeros(nRlz, 1);

% This stores the vector SupG for each run
supG_raw_store        = zeros(nBoot, nRlz);

%-These matrices store all the sets of interest during the bootstrap
% method for all levels of smoothing
lower_contour_raw_80_store       = zeros([nRlz dim]);
upper_contour_raw_80_store       = zeros([nRlz dim]);
upper_subset_mid_raw_80_store = zeros([nRlz dim]);
mid_subset_lower_raw_80_store = zeros([nRlz dim]);

lower_contour_raw_90_store       = zeros([nRlz dim]);
upper_contour_raw_90_store       = zeros([nRlz dim]);
upper_subset_mid_raw_90_store = zeros([nRlz dim]);
mid_subset_lower_raw_90_store = zeros([nRlz dim]);

lower_contour_raw_95_store       = zeros([nRlz dim]);
upper_contour_raw_95_store       = zeros([nRlz dim]);
upper_subset_mid_raw_95_store = zeros([nRlz dim]);
mid_subset_lower_raw_95_store = zeros([nRlz dim]);

supG_raw         = zeros(nBoot,1);

raw_field_boundary_store        = zeros(dim(2), nBoot);

% Creating linearly increasing signal across columns
Sig = repmat(linspace(1, 3), dim(2), 1);

% Uncomment to look at the Signal
%imagesc(Sig); axis image; colorbar
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
         biasmystd(reshape(observed_data,[prod(dim) nSubj]),stdblk),...
           dim); 
       
      % Residuals
      resid = bsxfun(@minus,observed_data,observed_mean);
      resid = spdiags(1./reshape(observed_std, [prod(dim) 1]), 0,prod(dim),prod(dim))*reshape(resid,[prod(dim) nSubj]); 
      
      % Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          gaussian_rand        = normrnd(0,1,[nSubj, 1]);
          resid_bootstrap      = resid*spdiags(gaussian_rand, 0, nSubj, nSubj);
          resid_bootstrap      = reshape(resid_bootstrap, [dim nSubj]);
          resid_field          = sum(resid_bootstrap, 3)/sqrt(nSubj); 

          supG_raw(k)          = max(abs(resid_field(true_boundary)));
          
          if t==nRlz
              raw_field_boundary_store(:,k)         = resid_field(true_boundary);
          end 
      end
    middle_contour                = AC;
      
    supGa_raw_80                     = prctile(supG_raw, 80);
    supGa_raw_90                     = prctile(supG_raw, 90);
    supGa_raw_95                     = prctile(supG_raw, 95);
       
    lower_contour_raw_80             = observed_mean >= thr - supGa_raw_80*tau*observed_std;
    upper_contour_raw_80             = observed_mean >= thr + supGa_raw_80*tau*observed_std;
    mid_on_upper_raw_80              = upper_contour_raw_80.*middle_contour;
    lower_on_mid_raw_80              = middle_contour.*lower_contour_raw_80;
    upper_subset_mid_raw_80          = upper_contour_raw_80 - mid_on_upper_raw_80;
    mid_subset_lower_raw_80          = middle_contour - lower_on_mid_raw_80;
    
    lower_contour_raw_80_boundary    = observed_mean(true_boundary) - supGa_raw_80*tau*observed_std(true_boundary);
    upper_contour_raw_80_boundary    = observed_mean(true_boundary) + supGa_raw_80*tau*observed_std(true_boundary);
    mid_on_upper_raw_80_boundary     = thr > upper_contour_raw_80_boundary;
    lower_on_mid_raw_80_boundary     = lower_contour_raw_80_boundary > thr;
    
    
    lower_contour_raw_90             = observed_mean >= thr - supGa_raw_90*tau*observed_std;
    upper_contour_raw_90             = observed_mean >= thr + supGa_raw_90*tau*observed_std;
    mid_on_upper_raw_90              = upper_contour_raw_90.*middle_contour;
    lower_on_mid_raw_90              = middle_contour.*lower_contour_raw_90;
    upper_subset_mid_raw_90          = upper_contour_raw_90 - mid_on_upper_raw_90;
    mid_subset_lower_raw_90          = middle_contour - lower_on_mid_raw_90;
    
    lower_contour_raw_90_boundary    = observed_mean(true_boundary) - supGa_raw_90*tau*observed_std(true_boundary);
    upper_contour_raw_90_boundary    = observed_mean(true_boundary) + supGa_raw_90*tau*observed_std(true_boundary);
    mid_on_upper_raw_90_boundary     = thr > upper_contour_raw_90_boundary;
    lower_on_mid_raw_90_boundary     = lower_contour_raw_90_boundary > thr;
    
    lower_contour_raw_95             = observed_mean >= thr - supGa_raw_95*tau*observed_std;
    upper_contour_raw_95             = observed_mean >= thr + supGa_raw_95*tau*observed_std;
    mid_on_upper_raw_95              = upper_contour_raw_95.*middle_contour;
    lower_on_mid_raw_95              = middle_contour.*lower_contour_raw_95;
    upper_subset_mid_raw_95          = upper_contour_raw_95 - mid_on_upper_raw_95;
    mid_subset_lower_raw_95          = middle_contour - lower_on_mid_raw_95; 
    
    lower_contour_raw_95_boundary    = observed_mean(true_boundary) - supGa_raw_95*tau*observed_std(true_boundary);
    upper_contour_raw_95_boundary    = observed_mean(true_boundary) + supGa_raw_95*tau*observed_std(true_boundary);
    mid_on_upper_raw_95_boundary     = thr > upper_contour_raw_95_boundary;
    lower_on_mid_raw_95_boundary     = lower_contour_raw_95_boundary > thr;
    
    %
    % Storing all variables of interest
    %
    supG_raw_store(:,t)                 = supG_raw;
    
    threshold_raw_80_store(t)              = supGa_raw_80;
    lower_contour_raw_80_store(t,:,:)      = lower_contour_raw_80;
    upper_contour_raw_80_store(t,:,:)      = upper_contour_raw_80;
    upper_subset_mid_raw_80_store(t,:,:)   = upper_subset_mid_raw_80;
    mid_subset_lower_raw_80_store(t,:,:)   = mid_subset_lower_raw_80;
    
    threshold_raw_90_store(t)              = supGa_raw_90;
    lower_contour_raw_90_store(t,:,:)      = lower_contour_raw_90;
    upper_contour_raw_90_store(t,:,:)      = upper_contour_raw_90;
    upper_subset_mid_raw_90_store(t,:,:)   = upper_subset_mid_raw_90;
    mid_subset_lower_raw_90_store(t,:,:)   = mid_subset_lower_raw_90;
    
    threshold_raw_95_store(t)              = supGa_raw_95;
    lower_contour_raw_95_store(t,:,:)      = lower_contour_raw_95;
    upper_contour_raw_95_store(t,:,:)      = upper_contour_raw_95;
    upper_subset_mid_raw_95_store(t,:,:)   = upper_subset_mid_raw_95;
    mid_subset_lower_raw_95_store(t,:,:)   = mid_subset_lower_raw_95;
    
    if sum(upper_subset_mid_raw_80(:))+sum(mid_subset_lower_raw_80(:))==0
      subset_success_vector_raw_80(t) = 1; 
      fprintf('raw nominal 80 success! \n');
    else 
      subset_success_vector_raw_80(t) = 0; 
      fprintf('raw nominal 80 failure! \n');
    end 
    
    if sum(mid_on_upper_raw_80_boundary(:))+sum(lower_on_mid_raw_80_boundary(:))==0
      subset_success_vector_raw_80_boundary(t) = 1;
      fprintf('raw nominal 80 boundary success! \n');
    else 
      subset_success_vector_raw_80_boundary(t) = 0;
      fprintf('raw nominal 80 boundary failure! \n');
    end 
    
    if sum(upper_subset_mid_raw_90(:))+sum(mid_subset_lower_raw_90(:))==0
      subset_success_vector_raw_90(t) = 1; 
      fprintf('raw nominal 90 success! \n');
    else 
      subset_success_vector_raw_90(t) = 0; 
      fprintf('raw nominal 90 failure! \n');
    end 
    
    if sum(mid_on_upper_raw_90_boundary(:))+sum(lower_on_mid_raw_90_boundary(:))==0
      subset_success_vector_raw_90_boundary(t) = 1;
      fprintf('raw nominal 90 boundary success! \n');
    else 
      subset_success_vector_raw_90_boundary(t) = 0;
      fprintf('raw nominal 90 boundary failure! \n');
    end 
    
    if sum(upper_subset_mid_raw_95(:))+sum(mid_subset_lower_raw_95(:))==0
      subset_success_vector_raw_95(t) = 1; 
      fprintf('raw nominal 95 success! \n');
    else 
      subset_success_vector_raw_95(t) = 0; 
      fprintf('raw nominal 95 failure! \n');
    end 
  
    if sum(mid_on_upper_raw_95_boundary(:))+sum(lower_on_mid_raw_95_boundary(:))==0
      subset_success_vector_raw_95_boundary(t) = 1;
      fprintf('raw nominal 95 boundary success! \n');
    else 
      subset_success_vector_raw_95_boundary(t) = 0;
      fprintf('raw nominal 95 boundary failure! \n');
    end 
end

percentage_success_vector_raw_80         = mean(subset_success_vector_raw_80, 1);
percentage_success_vector_raw_90         = mean(subset_success_vector_raw_90, 1);
percentage_success_vector_raw_95         = mean(subset_success_vector_raw_95, 1);
percentage_success_vector_raw_80_boundary         = mean(subset_success_vector_raw_80_boundary, 1);
percentage_success_vector_raw_90_boundary         = mean(subset_success_vector_raw_90_boundary, 1);
percentage_success_vector_raw_95_boundary         = mean(subset_success_vector_raw_95_boundary, 1);

eval(['save ' SvNm ' nSubj nRlz dim smo mag rimFWHM thr nBoot threshold_raw_80_store threshold_raw_90_store threshold_raw_95_store lower_contour_raw_80_store lower_contour_raw_90_store lower_contour_raw_95_store upper_contour_raw_80_store upper_contour_raw_90_store upper_contour_raw_95_store upper_subset_mid_raw_80_store upper_subset_mid_raw_90_store upper_subset_mid_raw_95_store mid_subset_lower_raw_80_store mid_subset_lower_raw_90_store mid_subset_lower_raw_95_store subset_success_vector_raw_80 subset_success_vector_raw_80_boundary subset_success_vector_raw_90 subset_success_vector_raw_90_boundary subset_success_vector_raw_95 subset_success_vector_raw_95_boundary percentage_success_vector_raw_80 percentage_success_vector_raw_80_boundary percentage_success_vector_raw_90 percentage_success_vector_raw_90_boundary percentage_success_vector_raw_95 percentage_success_vector_raw_95_boundary supG_raw_store raw_field_boundary_store '])
