function Sim_20(nSubj,SvNm,nRlz)
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
rad     = 30;

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


%- This vector stores the threshold value 'c' for each run
threshold_raw_80_store                  = zeros(nRlz, 1);
threshold_raw_90_store                  = zeros(nRlz, 1);
threshold_raw_95_store                  = zeros(nRlz, 1);

%- This vector stores the percentage volumes A^+_c/A_c, A^_c/A_c, A^-_c/A_c
lower_contour_raw_80_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_80_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_raw_90_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_90_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_raw_95_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_raw_95_volume_prct_store                     = zeros(nRlz, 1);

% This stores the vector SupG for each run
supG_raw_store                   = zeros(nBoot, nRlz);

%-These matrices store all the sets of interest during the bootstrap
% method for all levels of smoothing
lower_contour_raw_80_store                       = zeros([nRlz dim]);
upper_contour_raw_80_store                       = zeros([nRlz dim]);
upper_subset_mid_raw_80_store                    = zeros([nRlz dim]);
mid_subset_lower_raw_80_store                    = zeros([nRlz dim]);

lower_contour_raw_90_store                       = zeros([nRlz dim]);
upper_contour_raw_90_store                       = zeros([nRlz dim]);
upper_subset_mid_raw_90_store                    = zeros([nRlz dim]);
mid_subset_lower_raw_90_store                    = zeros([nRlz dim]);

lower_contour_raw_95_store                       = zeros([nRlz dim]);
upper_contour_raw_95_store                       = zeros([nRlz dim]);
upper_subset_mid_raw_95_store                    = zeros([nRlz dim]);
mid_subset_lower_raw_95_store                    = zeros([nRlz dim]);

supG_raw                         = zeros(nBoot,1);

% Creating linearly increasing signal across columns
Sig = repmat(linspace(1, 3), dim(2), 1);

% Uncomment to look at the Signal
%imagesc(Sig); axis image; colorbar
AC = Sig >= thr;

% Variables for computing the estimated boundary
[a,b] = ndgrid(-1:1);
se = strel('arbitrary',sqrt(a.^2 + b.^2) <=1);

% Obtaining the edges for the boundary Sig > 2 using the linear interpolation methods 
  % Making the interpolated boundary edges
  % Horizontal edges
  horz = AC(:,2:end) | AC(:,1:end-1);
  % Compute the left shifted horizontal edges
  lshift            = AC; % initialize
  lshift(:,1:end-1) = horz;
  lshift            = lshift & ~AC;
  %%% Compute the right shifted horizontal edges
  rshift          = AC; % initialize
  rshift(:,2:end) = horz;
  rshift          = rshift & ~AC;
  % Vertical edges
  vert = AC(1:end-1,:) | AC(2:end,:);
  %%% Compute the right shifted horizontal edges
  ushift = AC;
  ushift(1:end-1,:) = vert;
  ushift = ushift & ~AC;
  %%% Compute the down shifted vertical edges
  dshift = AC;
  %%% Values of random field on down shifted vertical edges
  dshift(2:end,:)   = vert;
  dshift = dshift & ~AC;
  
% Computing the weights for the weighted linear boundary
lshift_w1 = abs(Sig(lshift(:,[dim(2) 1:dim(2)-1])) - thr)./abs(Sig(lshift) - Sig(lshift(:,[dim(2) 1:dim(2)-1])));
lshift_w2 = abs(Sig(lshift) - thr)./abs(Sig(lshift) - Sig(lshift(:,[dim(2) 1:dim(2)-1])));

rshift_w1 = abs(Sig(rshift(:,[2:dim(2) 1])) - thr)./abs(Sig(rshift) - Sig(rshift(:,[2:dim(2) 1])));
rshift_w2 = abs(Sig(rshift) - thr)./abs(Sig(rshift) - Sig(rshift(:,[2:dim(2) 1])));

ushift_w1 = abs(Sig(ushift([dim(1) 1:dim(1)-1],:)) - thr)./abs(Sig(ushift) - Sig(ushift([dim(1) 1:dim(1)-1],:)));
ushift_w2 = abs(Sig(ushift) - thr)./abs(Sig(ushift) - Sig(ushift([dim(1) 1:dim(1)-1],:)));

dshift_w1 = abs(Sig(dshift([2:dim(1) 1],:)) - thr)./abs(Sig(dshift) - Sig(dshift([2:dim(1) 1],:)));
dshift_w2 = abs(Sig(dshift) - thr)./abs(Sig(dshift) - Sig(dshift([2:dim(1) 1],:)));


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
          % Applying the bootstrap using Rademacher variables (signflips)
          signflips                              = randi(2,[nSubj,1])*2-3;
          resid_bootstrap                        = resid*spdiags(signflips, 0, nSubj, nSubj);
          resid_bootstrap                        = reshape(resid_bootstrap, [dim nSubj]);
          resid_field                            = sum(resid_bootstrap, 3)/sqrt(nSubj); 

          % Calculating the maximum over the linear true weighted boundary edges
          lshift_boundary_values = abs(lshift_w1.*resid_field(lshift) + lshift_w2.*resid_field(lshift(:,[dim(2) 1:dim(2)-1])));
          rshift_boundary_values = abs(rshift_w1.*resid_field(rshift) + rshift_w2.*resid_field(rshift(:,[2:dim(2) 1])));
          ushift_boundary_values = abs(ushift_w1.*resid_field(ushift) + ushift_w2.*resid_field(ushift([dim(1) 1:dim(1)-1],:)));
          dshift_boundary_values = abs(dshift_w1.*resid_field(dshift) + dshift_w2.*resid_field(dshift([2:dim(1) 1],:)));
          supG_raw(k)            = max([lshift_boundary_values; rshift_boundary_values; ushift_boundary_values; dshift_boundary_values]);
          
      end
      
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

    
    lshift_observed_mean =  lshift_w1.*observed_mean(lshift) + lshift_w2.*observed_mean(lshift(:,[dim(2) 1:dim(2)-1]))

end

percentage_success_vector_raw_80                         = mean(subset_success_vector_raw_80, 1);
percentage_success_vector_raw_90                         = mean(subset_success_vector_raw_90, 1);
percentage_success_vector_raw_95                         = mean(subset_success_vector_raw_95, 1);


eval(['save ' SvNm ' nSubj nRlz dim smo mag rimFWHM thr nBoot '... 
      'threshold_raw_80_store threshold_raw_90_store threshold_raw_95_store threshold_raw_80_weighted_store threshold_raw_90_weighted_store threshold_raw_95_weighted_store threshold_raw_80_linear_store threshold_raw_90_linear_store threshold_raw_95_linear_store threshold_raw_80_observed_weighted_store threshold_raw_90_observed_weighted_store threshold_raw_95_observed_weighted_store '...
      'lower_contour_raw_80_store lower_contour_raw_90_store lower_contour_raw_95_store lower_contour_raw_80_weighted_store lower_contour_raw_90_weighted_store lower_contour_raw_95_weighted_store lower_contour_raw_80_linear_store lower_contour_raw_90_linear_store lower_contour_raw_95_linear_store lower_contour_raw_80_observed_weighted_store lower_contour_raw_90_observed_weighted_store lower_contour_raw_95_observed_weighted_store '...
      'upper_contour_raw_80_store upper_contour_raw_90_store upper_contour_raw_95_store upper_contour_raw_80_weighted_store upper_contour_raw_90_weighted_store upper_contour_raw_95_weighted_store upper_contour_raw_80_linear_store upper_contour_raw_90_linear_store upper_contour_raw_95_linear_store upper_contour_raw_80_observed_weighted_store upper_contour_raw_90_observed_weighted_store upper_contour_raw_95_observed_weighted_store '...
      'upper_subset_mid_raw_80_store upper_subset_mid_raw_90_store upper_subset_mid_raw_95_store upper_subset_mid_raw_80_weighted_store upper_subset_mid_raw_90_weighted_store upper_subset_mid_raw_95_weighted_store upper_subset_mid_raw_80_linear_store upper_subset_mid_raw_90_linear_store upper_subset_mid_raw_95_linear_store upper_subset_mid_raw_80_observed_weighted_store upper_subset_mid_raw_90_observed_weighted_store upper_subset_mid_raw_95_observed_weighted_store '...
      'mid_subset_lower_raw_80_store mid_subset_lower_raw_90_store mid_subset_lower_raw_95_store mid_subset_lower_raw_80_weighted_store mid_subset_lower_raw_90_weighted_store mid_subset_lower_raw_95_weighted_store mid_subset_lower_raw_80_linear_store mid_subset_lower_raw_90_linear_store mid_subset_lower_raw_95_linear_store mid_subset_lower_raw_80_observed_weighted_store mid_subset_lower_raw_90_observed_weighted_store mid_subset_lower_raw_95_observed_weighted_store '...
      'subset_success_vector_raw_80 subset_success_vector_raw_90 subset_success_vector_raw_95 subset_success_vector_raw_80_weighted subset_success_vector_raw_90_weighted subset_success_vector_raw_95_weighted subset_success_vector_raw_80_linear subset_success_vector_raw_90_linear subset_success_vector_raw_95_linear subset_success_vector_raw_80_observed_weighted subset_success_vector_raw_90_observed_weighted subset_success_vector_raw_95_observed_weighted '...
      'percentage_success_vector_raw_80 percentage_success_vector_raw_90 percentage_success_vector_raw_95 percentage_success_vector_raw_80_weighted percentage_success_vector_raw_90_weighted percentage_success_vector_raw_95_weighted percentage_success_vector_raw_80_linear percentage_success_vector_raw_90_linear percentage_success_vector_raw_95_linear percentage_success_vector_raw_80_observed_weighted percentage_success_vector_raw_90_observed_weighted percentage_success_vector_raw_95_observed_weighted '...
      'supG_raw_store supG_raw_weighted_store supG_raw_linear_store supG_raw_observed_weighted_store '...
      'middle_contour_volume observed_AC_volume '...
      'lower_contour_raw_80_volume_prct_store lower_contour_raw_90_volume_prct_store lower_contour_raw_95_volume_prct_store lower_contour_raw_80_weighted_volume_prct_store lower_contour_raw_90_weighted_volume_prct_store lower_contour_raw_95_weighted_volume_prct_store lower_contour_raw_80_linear_volume_prct_store lower_contour_raw_90_linear_volume_prct_store lower_contour_raw_95_linear_volume_prct_store lower_contour_raw_80_observed_weighted_volume_prct_store lower_contour_raw_90_observed_weighted_volume_prct_store lower_contour_raw_95_observed_weighted_volume_prct_store '...
      'upper_contour_raw_80_volume_prct_store upper_contour_raw_90_volume_prct_store upper_contour_raw_95_volume_prct_store upper_contour_raw_80_weighted_volume_prct_store upper_contour_raw_90_weighted_volume_prct_store upper_contour_raw_95_weighted_volume_prct_store upper_contour_raw_80_linear_volume_prct_store upper_contour_raw_90_linear_volume_prct_store upper_contour_raw_95_linear_volume_prct_store upper_contour_raw_80_observed_weighted_volume_prct_store upper_contour_raw_90_observed_weighted_volume_prct_store upper_contour_raw_95_observed_weighted_volume_prct_store'])
