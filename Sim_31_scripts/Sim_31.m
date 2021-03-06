function Sim_31(nSubj,SvNm,nRlz)
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
subset_success_vector_observed_80           = zeros(nRlz, 1); 
subset_success_vector_observed_90           = zeros(nRlz, 1);
subset_success_vector_observed_95           = zeros(nRlz, 1);
subset_success_vector_cohen_d_80           = zeros(nRlz, 1); 
subset_success_vector_cohen_d_90           = zeros(nRlz, 1);
subset_success_vector_cohen_d_95           = zeros(nRlz, 1);
subset_success_vector_observed_80_alternate = zeros(nRlz, 1); 
subset_success_vector_observed_90_alternate = zeros(nRlz, 1);
subset_success_vector_observed_95_alternate = zeros(nRlz, 1);
subset_success_vector_cohen_d_80_alternate = zeros(nRlz, 1); 
subset_success_vector_cohen_d_90_alternate = zeros(nRlz, 1);
subset_success_vector_cohen_d_95_alternate = zeros(nRlz, 1);

%- This vector stores the threshold value 'c' for each run
threshold_observed_80_store                  = zeros(nRlz, 1);
threshold_observed_90_store                  = zeros(nRlz, 1);
threshold_observed_95_store                  = zeros(nRlz, 1);

threshold_cohen_d_80_store                  = zeros(nRlz, 1);
threshold_cohen_d_90_store                  = zeros(nRlz, 1);
threshold_cohen_d_95_store                  = zeros(nRlz, 1);

%- This vector stores the percentage volumes A^+_c/A_c, A^_c/A_c, A^-_c/A_c
lower_contour_observed_80_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_80_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_observed_90_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_90_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_observed_95_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_observed_95_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_cohen_d_80_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_cohen_d_80_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_cohen_d_90_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_cohen_d_90_volume_prct_store                     = zeros(nRlz, 1);

lower_contour_cohen_d_95_volume_prct_store                     = zeros(nRlz, 1);
upper_contour_cohen_d_95_volume_prct_store                     = zeros(nRlz, 1);

% This stores the vector SupG for each run
supG_observed_store              = zeros(nBoot, nRlz);
supG_cohen_d_store               = zeros(nBoot, nRlz);

%-These matrices store all the sets of interest during the bootstrap
% method for all levels of smoothing
lower_contour_observed_80_store                       = zeros([nRlz dim]);
upper_contour_observed_80_store                       = zeros([nRlz dim]);
upper_subset_mid_observed_80_store                    = zeros([nRlz dim]);
mid_subset_lower_observed_80_store                    = zeros([nRlz dim]);

lower_contour_observed_90_store                       = zeros([nRlz dim]);
upper_contour_observed_90_store                       = zeros([nRlz dim]);
upper_subset_mid_observed_90_store                    = zeros([nRlz dim]);
mid_subset_lower_observed_90_store                    = zeros([nRlz dim]);

lower_contour_observed_95_store                       = zeros([nRlz dim]);
upper_contour_observed_95_store                       = zeros([nRlz dim]);
upper_subset_mid_observed_95_store                    = zeros([nRlz dim]);
mid_subset_lower_observed_95_store                    = zeros([nRlz dim]);

lower_contour_cohen_d_80_store                       = zeros([nRlz dim]);
upper_contour_cohen_d_80_store                       = zeros([nRlz dim]);
upper_subset_mid_cohen_d_80_store                    = zeros([nRlz dim]);
mid_subset_lower_cohen_d_80_store                    = zeros([nRlz dim]);

lower_contour_cohen_d_90_store                       = zeros([nRlz dim]);
upper_contour_cohen_d_90_store                       = zeros([nRlz dim]);
upper_subset_mid_cohen_d_90_store                    = zeros([nRlz dim]);
mid_subset_lower_cohen_d_90_store                    = zeros([nRlz dim]);

lower_contour_cohen_d_95_store                       = zeros([nRlz dim]);
upper_contour_cohen_d_95_store                       = zeros([nRlz dim]);
upper_subset_mid_cohen_d_95_store                    = zeros([nRlz dim]);
mid_subset_lower_cohen_d_95_store                    = zeros([nRlz dim]);

supG_observed                    = zeros(nBoot,1);
supG_cohen_d                     = zeros(nBoot,1);

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
      
      % Making the observed boundary and weights
      observed_AC = observed_mean >= thr;
      
      % Making the interpolated boundary edges
      % Horizontal edges
      observed_horz = observed_AC(:,2:end) | observed_AC(:,1:end-1);
      % Compute the left shifted horizontal edges
      observed_lshift            = observed_AC; % initialize
      observed_lshift(:,1:end-1) = observed_horz;
      observed_lshift            = observed_lshift & ~observed_AC;
      %%% Compute the right shifted horizontal edges
      observed_rshift          = observed_AC; % initialize
      observed_rshift(:,2:end) = observed_horz;
      observed_rshift          = observed_rshift & ~observed_AC;
      % Vertical edges
      vert = observed_AC(1:end-1,:) | observed_AC(2:end,:);
      %%% Compute the up shifted horizontal edges
      observed_ushift = observed_AC;
      observed_ushift(1:end-1,:) = vert;
      observed_ushift = observed_ushift & ~observed_AC;
      %%% Compute the down shifted vertical edges
      observed_dshift = observed_AC;
      observed_dshift(2:end,:)   = vert;
      observed_dshift = observed_dshift & ~observed_AC;
      
      % Computing the weights for the weighted linear boundary
      observed_lshift_w1 = abs(observed_mean(observed_lshift(:,[dim(2) 1:dim(2)-1])) - thr)./abs(observed_mean(observed_lshift) - observed_mean(observed_lshift(:,[dim(2) 1:dim(2)-1])));
      observed_lshift_w2 = abs(observed_mean(observed_lshift) - thr)./abs(observed_mean(observed_lshift) - observed_mean(observed_lshift(:,[dim(2) 1:dim(2)-1])));

      observed_rshift_w1 = abs(observed_mean(observed_rshift(:,[2:dim(2) 1])) - thr)./abs(observed_mean(observed_rshift) - observed_mean(observed_rshift(:,[2:dim(2) 1])));
      observed_rshift_w2 = abs(observed_mean(observed_rshift) - thr)./abs(observed_mean(observed_rshift) - observed_mean(observed_rshift(:,[2:dim(2) 1])));

      observed_ushift_w1 = abs(observed_mean(observed_ushift([dim(1) 1:dim(1)-1],:)) - thr)./abs(observed_mean(observed_ushift) - observed_mean(observed_ushift([dim(1) 1:dim(1)-1],:)));
      observed_ushift_w2 = abs(observed_mean(observed_ushift) - thr)./abs(observed_mean(observed_ushift) - observed_mean(observed_ushift([dim(1) 1:dim(1)-1],:)));

      observed_dshift_w1 = abs(observed_mean(observed_dshift([2:dim(1) 1],:)) - thr)./abs(observed_mean(observed_dshift) - observed_mean(observed_dshift([2:dim(1) 1],:)));
      observed_dshift_w2 = abs(observed_mean(observed_dshift) - thr)./abs(observed_mean(observed_dshift) - observed_mean(observed_dshift([2:dim(1) 1],:)));

      % Residuals
      resid = bsxfun(@minus,observed_data,observed_mean);
      resid = spdiags(1./reshape(observed_std, [prod(dim) 1]), 0,prod(dim),prod(dim))*reshape(resid,[prod(dim) nSubj]);
            
      %% Computing variables of interest for the Cohen's d case
      cohen_d = observed_mean./observed_std;
      cohen_d_std      = sqrt(1 + cohen_d.^2/2);

      cohen_d_AC       = cohen_d >= thr;
      
      % Calculating boundary edges for cohen d
      % Making the interpolated boundary edges
      % Horizontal edges
      cohen_d_horz                  = cohen_d_AC(:,2:end,:) | cohen_d_AC(:,1:end-1,:);
      % Compute the left shifted horizontal edges
      cohen_d_lshift               = cohen_d_AC; % initialize
      cohen_d_lshift(:,1:end-1,:)  = cohen_d_horz;
      cohen_d_lshift               = cohen_d_lshift & ~cohen_d_AC;
      %%% Compute the right shifted horizontal edges
      cohen_d_rshift               = cohen_d_AC; % initialize
      cohen_d_rshift(:,2:end,:)    = cohen_d_horz;
      cohen_d_rshift               = cohen_d_rshift & ~cohen_d_AC;
      % Vertical edges
      cohen_d_vert			      = cohen_d_AC(1:end-1,:,:) | cohen_d_AC(2:end,:,:);
      %%% Compute the up shifted horizontal edges
      cohen_d_ushift               = cohen_d_AC;
      cohen_d_ushift(1:end-1,:,:)  = cohen_d_vert;
      cohen_d_ushift               = cohen_d_ushift & ~cohen_d_AC;
      %%% Compute the down shifted vertical edges
      cohen_d_dshift               = cohen_d_AC;
      cohen_d_dshift(2:end,:,:)    = cohen_d_vert;
      cohen_d_dshift               = cohen_d_dshift & ~cohen_d_AC;

      % Computing the weights for the weighted linear boundary
      cohen_d_lshift_w1 = abs(cohen_d(cohen_d_lshift(:,[dim(2) 1:dim(2)-1],:)) - thr)./abs(cohen_d(cohen_d_lshift) - cohen_d(cohen_d_lshift(:,[dim(2) 1:dim(2)-1],:)));
      cohen_d_lshift_w2 = abs(cohen_d(cohen_d_lshift) - thr)./abs(cohen_d(cohen_d_lshift) - cohen_d(cohen_d_lshift(:,[dim(2) 1:dim(2)-1],:)));

      cohen_d_rshift_w1 = abs(cohen_d(cohen_d_rshift(:,[2:dim(2) 1],:)) - thr)./abs(cohen_d(cohen_d_rshift) - cohen_d(cohen_d_rshift(:,[2:dim(2) 1],:)));
      cohen_d_rshift_w2 = abs(cohen_d(cohen_d_rshift) - thr)./abs(cohen_d(cohen_d_rshift) - cohen_d(cohen_d_rshift(:,[2:dim(2) 1],:)));

      cohen_d_ushift_w1 = abs(cohen_d(cohen_d_ushift([dim(1) 1:dim(1)-1],:,:)) - thr)./abs(cohen_d(cohen_d_ushift) - cohen_d(cohen_d_ushift([dim(1) 1:dim(1)-1],:,:)));
      cohen_d_ushift_w2 = abs(cohen_d(cohen_d_ushift) - thr)./abs(cohen_d(cohen_d_ushift) - cohen_d(cohen_d_ushift([dim(1) 1:dim(1)-1],:,:)));

      cohen_d_dshift_w1 = abs(cohen_d(cohen_d_dshift([2:dim(1) 1],:,:)) - thr)./abs(cohen_d(cohen_d_dshift) - cohen_d(cohen_d_dshift([2:dim(1) 1],:,:)));
      cohen_d_dshift_w2 = abs(cohen_d(cohen_d_dshift) - thr)./abs(cohen_d(cohen_d_dshift) - cohen_d(cohen_d_dshift([2:dim(1) 1],:,:)));
      
      % Residuals using the SNR transformation
      cohen_d_resid = ( (bsxfun(@minus, observed_data, observed_mean))./observed_std - cohen_d/2.*(((bsxfun(@minus, observed_data, observed_mean))./observed_std).^2-1) )...
                    ./ cohen_d_std;
      cohen_d_resid = reshape(cohen_d_resid, [prod(dim), nSubj]);
      
      % Implementing the Multiplier Boostrap to obtain confidence intervals
      for k=1:nBoot 
          % Applying the bootstrap using Rademacher variables (signflips)
          signflips                              = randi(2,[nSubj,1])*2-3;
          resid_bootstrap                        = resid*spdiags(signflips, 0, nSubj, nSubj);
          resid_bootstrap                        = reshape(resid_bootstrap, [dim nSubj]);
          resid_field                            = sum(resid_bootstrap, 3)/sqrt(nSubj); 

          % Calculating the maximum over the mu^ = thresh weighted boundary
          % edges
          observed_lshift_boundary_values = abs(observed_lshift_w1.*resid_field(observed_lshift) + observed_lshift_w2.*resid_field(observed_lshift(:,[dim(2) 1:dim(2)-1])));
          observed_rshift_boundary_values = abs(observed_rshift_w1.*resid_field(observed_rshift) + observed_rshift_w2.*resid_field(observed_rshift(:,[2:dim(2) 1])));
          observed_ushift_boundary_values = abs(observed_ushift_w1.*resid_field(observed_ushift) + observed_ushift_w2.*resid_field(observed_ushift([dim(1) 1:dim(1)-1],:)));
          observed_dshift_boundary_values = abs(observed_dshift_w1.*resid_field(observed_dshift) + observed_dshift_w2.*resid_field(observed_dshift([2:dim(1) 1],:)));
          supG_observed(k)            = max([observed_lshift_boundary_values; observed_rshift_boundary_values; observed_ushift_boundary_values; observed_dshift_boundary_values]);
          
          % Cohen's d case
          cohen_d_resid_bootstrap                        = cohen_d_resid*spdiags(signflips, 0, nSubj, nSubj);
          cohen_d_resid_bootstrap                        = reshape(cohen_d_resid_bootstrap, [dim nSubj]);
          cohen_d_resid_field                            = sum(cohen_d_resid_bootstrap, 3)/sqrt(nSubj); 
 
          % Calculating the maximum over the mu^ = thresh weighted boundary
          % edges
          cohen_d_lshift_boundary_values = abs(cohen_d_lshift_w1.*cohen_d_resid_field(cohen_d_lshift) + cohen_d_lshift_w2.*cohen_d_resid_field(cohen_d_lshift(:,[dim(2) 1:dim(2)-1])));
          cohen_d_rshift_boundary_values = abs(cohen_d_rshift_w1.*cohen_d_resid_field(cohen_d_rshift) + cohen_d_rshift_w2.*cohen_d_resid_field(cohen_d_rshift(:,[2:dim(2) 1])));
          cohen_d_ushift_boundary_values = abs(cohen_d_ushift_w1.*cohen_d_resid_field(cohen_d_ushift) + cohen_d_ushift_w2.*cohen_d_resid_field(cohen_d_ushift([dim(1) 1:dim(1)-1],:)));
          cohen_d_dshift_boundary_values = abs(cohen_d_dshift_w1.*cohen_d_resid_field(cohen_d_dshift) + cohen_d_dshift_w2.*cohen_d_resid_field(cohen_d_dshift([2:dim(1) 1],:)));
          supG_cohen_d(k)            = max([cohen_d_lshift_boundary_values; cohen_d_rshift_boundary_values; cohen_d_ushift_boundary_values; cohen_d_dshift_boundary_values]);
          
      end
      
    middle_contour                = AC;
    middle_contour_volume         = sum(middle_contour(:));
    
    % Gaussian random variable results for the true and estimated boundary
    % Observed boundary
    supGa_observed_80                     = prctile(supG_observed, 80);
    supGa_observed_90                     = prctile(supG_observed, 90);
    supGa_observed_95                     = prctile(supG_observed, 95);
       
    lower_contour_observed_80             = observed_mean >= thr - supGa_observed_80*tau*observed_std;
    upper_contour_observed_80             = observed_mean >= thr + supGa_observed_80*tau*observed_std;
    lower_contour_observed_80_volume_prct = sum(lower_contour_observed_80(:))/middle_contour_volume;
    upper_contour_observed_80_volume_prct = sum(upper_contour_observed_80(:))/middle_contour_volume;
    mid_on_upper_observed_80              = upper_contour_observed_80.*middle_contour;
    lower_on_mid_observed_80              = middle_contour.*lower_contour_observed_80;
    upper_subset_mid_observed_80          = upper_contour_observed_80 - mid_on_upper_observed_80;
    mid_subset_lower_observed_80          = middle_contour - lower_on_mid_observed_80;
    
    lower_contour_observed_90             = observed_mean >= thr - supGa_observed_90*tau*observed_std;
    upper_contour_observed_90             = observed_mean >= thr + supGa_observed_90*tau*observed_std;
    lower_contour_observed_90_volume_prct = sum(lower_contour_observed_90(:))/middle_contour_volume;
    upper_contour_observed_90_volume_prct = sum(upper_contour_observed_90(:))/middle_contour_volume;
    mid_on_upper_observed_90              = upper_contour_observed_90.*middle_contour;
    lower_on_mid_observed_90              = middle_contour.*lower_contour_observed_90;
    upper_subset_mid_observed_90          = upper_contour_observed_90 - mid_on_upper_observed_90;
    mid_subset_lower_observed_90          = middle_contour - lower_on_mid_observed_90;    
    
    lower_contour_observed_95             = observed_mean >= thr - supGa_observed_95*tau*observed_std;
    upper_contour_observed_95             = observed_mean >= thr + supGa_observed_95*tau*observed_std;
    lower_contour_observed_95_volume_prct = sum(lower_contour_observed_95(:))/middle_contour_volume;
    upper_contour_observed_95_volume_prct = sum(upper_contour_observed_95(:))/middle_contour_volume;
    mid_on_upper_observed_95              = upper_contour_observed_95.*middle_contour;
    lower_on_mid_observed_95              = middle_contour.*lower_contour_observed_95;
    upper_subset_mid_observed_95          = upper_contour_observed_95 - mid_on_upper_observed_95;
    mid_subset_lower_observed_95          = middle_contour - lower_on_mid_observed_95;

    % Cohen's d results
    supGa_cohen_d_80                     = prctile(supG_cohen_d, 80);
    supGa_cohen_d_90                     = prctile(supG_cohen_d, 90);
    supGa_cohen_d_95                     = prctile(supG_cohen_d, 95);
       
    lower_contour_cohen_d_80             = cohen_d >= thr - supGa_cohen_d_80*tau*cohen_d_std;
    upper_contour_cohen_d_80             = cohen_d >= thr + supGa_cohen_d_80*tau*cohen_d_std;
    lower_contour_cohen_d_80_volume_prct = sum(lower_contour_cohen_d_80(:))/middle_contour_volume;
    upper_contour_cohen_d_80_volume_prct = sum(upper_contour_cohen_d_80(:))/middle_contour_volume;
    mid_on_upper_cohen_d_80              = upper_contour_cohen_d_80.*middle_contour;
    lower_on_mid_cohen_d_80              = middle_contour.*lower_contour_cohen_d_80;
    upper_subset_mid_cohen_d_80          = upper_contour_cohen_d_80 - mid_on_upper_cohen_d_80;
    mid_subset_lower_cohen_d_80          = middle_contour - lower_on_mid_cohen_d_80;
    
    lower_contour_cohen_d_90             = cohen_d >= thr - supGa_cohen_d_90*tau*cohen_d_std;
    upper_contour_cohen_d_90             = cohen_d >= thr + supGa_cohen_d_90*tau*cohen_d_std;
    lower_contour_cohen_d_90_volume_prct = sum(lower_contour_cohen_d_90(:))/middle_contour_volume;
    upper_contour_cohen_d_90_volume_prct = sum(upper_contour_cohen_d_90(:))/middle_contour_volume;
    mid_on_upper_cohen_d_90              = upper_contour_cohen_d_90.*middle_contour;
    lower_on_mid_cohen_d_90              = middle_contour.*lower_contour_cohen_d_90;
    upper_subset_mid_cohen_d_90          = upper_contour_cohen_d_90 - mid_on_upper_cohen_d_90;
    mid_subset_lower_cohen_d_90          = middle_contour - lower_on_mid_cohen_d_90;    
    
    lower_contour_cohen_d_95             = cohen_d >= thr - supGa_cohen_d_95*tau*cohen_d_std;
    upper_contour_cohen_d_95             = cohen_d >= thr + supGa_cohen_d_95*tau*cohen_d_std;
    lower_contour_cohen_d_95_volume_prct = sum(lower_contour_cohen_d_95(:))/middle_contour_volume;
    upper_contour_cohen_d_95_volume_prct = sum(upper_contour_cohen_d_95(:))/middle_contour_volume;
    mid_on_upper_cohen_d_95              = upper_contour_cohen_d_95.*middle_contour;
    lower_on_mid_cohen_d_95              = middle_contour.*lower_contour_cohen_d_95;
    upper_subset_mid_cohen_d_95          = upper_contour_cohen_d_95 - mid_on_upper_cohen_d_95;
    mid_subset_lower_cohen_d_95          = middle_contour - lower_on_mid_cohen_d_95;

    %
    % Storing all variables of interest
    %
    % Observed boundary variables
    supG_observed_store(:,t)                                    = supG_observed;
    threshold_observed_80_store(t)                              = supGa_observed_80;
    lower_contour_observed_80_store(t,:,:)                      = lower_contour_observed_80;
    upper_contour_observed_80_store(t,:,:)                      = upper_contour_observed_80;
    upper_subset_mid_observed_80_store(t,:,:)                   = upper_subset_mid_observed_80;
    mid_subset_lower_observed_80_store(t,:,:)                   = mid_subset_lower_observed_80;
    lower_contour_observed_80_volume_prct_store(t)              = lower_contour_observed_80_volume_prct;
    upper_contour_observed_80_volume_prct_store(t)              = upper_contour_observed_80_volume_prct;
 
    threshold_observed_90_store(t)                              = supGa_observed_90;
    lower_contour_observed_90_store(t,:,:)                      = lower_contour_observed_90;
    upper_contour_observed_90_store(t,:,:)                      = upper_contour_observed_90;
    upper_subset_mid_observed_90_store(t,:,:)                   = upper_subset_mid_observed_90;
    mid_subset_lower_observed_90_store(t,:,:)                   = mid_subset_lower_observed_90;
    lower_contour_observed_90_volume_prct_store(t)              = lower_contour_observed_90_volume_prct;
    upper_contour_observed_90_volume_prct_store(t)              = upper_contour_observed_90_volume_prct;

    threshold_observed_95_store(t)                              = supGa_observed_95;
    lower_contour_observed_95_store(t,:,:)                      = lower_contour_observed_95;
    upper_contour_observed_95_store(t,:,:)                      = upper_contour_observed_95;
    upper_subset_mid_observed_95_store(t,:,:)                   = upper_subset_mid_observed_95;
    mid_subset_lower_observed_95_store(t,:,:)                   = mid_subset_lower_observed_95;
    lower_contour_observed_95_volume_prct_store(t)              = lower_contour_observed_95_volume_prct;
    upper_contour_observed_95_volume_prct_store(t)              = upper_contour_observed_95_volume_prct;
    
    % Cohen d variables
    supG_cohen_d_store(:,t)                                    = supG_cohen_d;
    threshold_cohen_d_80_store(t)                              = supGa_cohen_d_80;
    lower_contour_cohen_d_80_store(t,:,:)                      = lower_contour_cohen_d_80;
    upper_contour_cohen_d_80_store(t,:,:)                      = upper_contour_cohen_d_80;
    upper_subset_mid_cohen_d_80_store(t,:,:)                   = upper_subset_mid_cohen_d_80;
    mid_subset_lower_cohen_d_80_store(t,:,:)                   = mid_subset_lower_cohen_d_80;
    lower_contour_cohen_d_80_volume_prct_store(t)              = lower_contour_cohen_d_80_volume_prct;
    upper_contour_cohen_d_80_volume_prct_store(t)              = upper_contour_cohen_d_80_volume_prct;
 
    threshold_cohen_d_90_store(t)                              = supGa_cohen_d_90;
    lower_contour_cohen_d_90_store(t,:,:)                      = lower_contour_cohen_d_90;
    upper_contour_cohen_d_90_store(t,:,:)                      = upper_contour_cohen_d_90;
    upper_subset_mid_cohen_d_90_store(t,:,:)                   = upper_subset_mid_cohen_d_90;
    mid_subset_lower_cohen_d_90_store(t,:,:)                   = mid_subset_lower_cohen_d_90;
    lower_contour_cohen_d_90_volume_prct_store(t)              = lower_contour_cohen_d_90_volume_prct;
    upper_contour_cohen_d_90_volume_prct_store(t)              = upper_contour_cohen_d_90_volume_prct;

    threshold_cohen_d_95_store(t)                              = supGa_cohen_d_95;
    lower_contour_cohen_d_95_store(t,:,:)                      = lower_contour_cohen_d_95;
    upper_contour_cohen_d_95_store(t,:,:)                      = upper_contour_cohen_d_95;
    upper_subset_mid_cohen_d_95_store(t,:,:)                   = upper_subset_mid_cohen_d_95;
    mid_subset_lower_cohen_d_95_store(t,:,:)                   = mid_subset_lower_cohen_d_95;
    lower_contour_cohen_d_95_volume_prct_store(t)              = lower_contour_cohen_d_95_volume_prct;
    upper_contour_cohen_d_95_volume_prct_store(t)              = upper_contour_cohen_d_95_volume_prct;
    
    % Sampling the observed mean process along the true boundary u = c
    lshift_observed_mean_boundary = abs(lshift_w1.*observed_mean(lshift) + lshift_w2.*observed_mean(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_observed_mean_boundary = abs(rshift_w1.*observed_mean(rshift) + rshift_w2.*observed_mean(rshift(:,[2:dim(2) 1])));
    ushift_observed_mean_boundary = abs(ushift_w1.*observed_mean(ushift) + ushift_w2.*observed_mean(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_observed_mean_boundary = abs(dshift_w1.*observed_mean(dshift) + dshift_w2.*observed_mean(dshift([2:dim(1) 1],:)));
    
    % Calculating the subset condition when residuals in multiplier
    % bootstrap are taken along the observed boundary
    lower_condition_80_observed = thr - supGa_observed_80*tau*observed_std;
    upper_condition_80_observed = thr + supGa_observed_80*tau*observed_std;
    lower_condition_90_observed = thr - supGa_observed_90*tau*observed_std;
    upper_condition_90_observed = thr + supGa_observed_90*tau*observed_std;
    lower_condition_95_observed = thr - supGa_observed_95*tau*observed_std;
    upper_condition_95_observed = thr + supGa_observed_95*tau*observed_std;
    
    lshift_lower_condition_80_observed_boundary = abs(lshift_w1.*lower_condition_80_observed(lshift) + lshift_w2.*lower_condition_80_observed(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_lower_condition_80_observed_boundary = abs(rshift_w1.*lower_condition_80_observed(rshift) + rshift_w2.*lower_condition_80_observed(rshift(:,[2:dim(2) 1])));
    ushift_lower_condition_80_observed_boundary = abs(ushift_w1.*lower_condition_80_observed(ushift) + ushift_w2.*lower_condition_80_observed(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_lower_condition_80_observed_boundary = abs(dshift_w1.*lower_condition_80_observed(dshift) + dshift_w2.*lower_condition_80_observed(dshift([2:dim(1) 1],:)));
    
    lshift_upper_condition_80_observed_boundary = abs(lshift_w1.*upper_condition_80_observed(lshift) + lshift_w2.*upper_condition_80_observed(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_upper_condition_80_observed_boundary = abs(rshift_w1.*upper_condition_80_observed(rshift) + rshift_w2.*upper_condition_80_observed(rshift(:,[2:dim(2) 1])));
    ushift_upper_condition_80_observed_boundary = abs(ushift_w1.*upper_condition_80_observed(ushift) + ushift_w2.*upper_condition_80_observed(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_upper_condition_80_observed_boundary = abs(dshift_w1.*upper_condition_80_observed(dshift) + dshift_w2.*upper_condition_80_observed(dshift([2:dim(1) 1],:)));    
    
    lshift_lower_condition_90_observed_boundary = abs(lshift_w1.*lower_condition_90_observed(lshift) + lshift_w2.*lower_condition_90_observed(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_lower_condition_90_observed_boundary = abs(rshift_w1.*lower_condition_90_observed(rshift) + rshift_w2.*lower_condition_90_observed(rshift(:,[2:dim(2) 1])));
    ushift_lower_condition_90_observed_boundary = abs(ushift_w1.*lower_condition_90_observed(ushift) + ushift_w2.*lower_condition_90_observed(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_lower_condition_90_observed_boundary = abs(dshift_w1.*lower_condition_90_observed(dshift) + dshift_w2.*lower_condition_90_observed(dshift([2:dim(1) 1],:)));
    
    lshift_upper_condition_90_observed_boundary = abs(lshift_w1.*upper_condition_90_observed(lshift) + lshift_w2.*upper_condition_90_observed(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_upper_condition_90_observed_boundary = abs(rshift_w1.*upper_condition_90_observed(rshift) + rshift_w2.*upper_condition_90_observed(rshift(:,[2:dim(2) 1])));
    ushift_upper_condition_90_observed_boundary = abs(ushift_w1.*upper_condition_90_observed(ushift) + ushift_w2.*upper_condition_90_observed(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_upper_condition_90_observed_boundary = abs(dshift_w1.*upper_condition_90_observed(dshift) + dshift_w2.*upper_condition_90_observed(dshift([2:dim(1) 1],:)));
    
    lshift_lower_condition_95_observed_boundary = abs(lshift_w1.*lower_condition_95_observed(lshift) + lshift_w2.*lower_condition_95_observed(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_lower_condition_95_observed_boundary = abs(rshift_w1.*lower_condition_95_observed(rshift) + rshift_w2.*lower_condition_95_observed(rshift(:,[2:dim(2) 1])));
    ushift_lower_condition_95_observed_boundary = abs(ushift_w1.*lower_condition_95_observed(ushift) + ushift_w2.*lower_condition_95_observed(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_lower_condition_95_observed_boundary = abs(dshift_w1.*lower_condition_95_observed(dshift) + dshift_w2.*lower_condition_95_observed(dshift([2:dim(1) 1],:)));
    
    lshift_upper_condition_95_observed_boundary = abs(lshift_w1.*upper_condition_95_observed(lshift) + lshift_w2.*upper_condition_95_observed(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_upper_condition_95_observed_boundary = abs(rshift_w1.*upper_condition_95_observed(rshift) + rshift_w2.*upper_condition_95_observed(rshift(:,[2:dim(2) 1])));
    ushift_upper_condition_95_observed_boundary = abs(ushift_w1.*upper_condition_95_observed(ushift) + ushift_w2.*upper_condition_95_observed(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_upper_condition_95_observed_boundary = abs(dshift_w1.*upper_condition_95_observed(dshift) + dshift_w2.*upper_condition_95_observed(dshift([2:dim(1) 1],:)));
    
    lower_condition_80_observed_success = [lshift_observed_mean_boundary < lshift_lower_condition_80_observed_boundary, ... 
                                  rshift_observed_mean_boundary < rshift_lower_condition_80_observed_boundary, ...
                                  ushift_observed_mean_boundary < ushift_lower_condition_80_observed_boundary, ...
                                  dshift_observed_mean_boundary < dshift_lower_condition_80_observed_boundary];
    upper_condition_80_observed_success = [lshift_observed_mean_boundary >= lshift_upper_condition_80_observed_boundary, ... 
                                  rshift_observed_mean_boundary >= rshift_upper_condition_80_observed_boundary, ...
                                  ushift_observed_mean_boundary >= ushift_upper_condition_80_observed_boundary, ...
                                  dshift_observed_mean_boundary >= dshift_upper_condition_80_observed_boundary];
                              
    lower_condition_90_observed_success = [lshift_observed_mean_boundary < lshift_lower_condition_90_observed_boundary, ... 
                                  rshift_observed_mean_boundary < rshift_lower_condition_90_observed_boundary, ...
                                  ushift_observed_mean_boundary < ushift_lower_condition_90_observed_boundary, ...
                                  dshift_observed_mean_boundary < dshift_lower_condition_90_observed_boundary];
    upper_condition_90_observed_success = [lshift_observed_mean_boundary >= lshift_upper_condition_90_observed_boundary, ... 
                                  rshift_observed_mean_boundary >= rshift_upper_condition_90_observed_boundary, ...
                                  ushift_observed_mean_boundary >= ushift_upper_condition_90_observed_boundary, ...
                                  dshift_observed_mean_boundary >= dshift_upper_condition_90_observed_boundary];
                              
    lower_condition_95_observed_success = [lshift_observed_mean_boundary < lshift_lower_condition_95_observed_boundary, ... 
                                  rshift_observed_mean_boundary < rshift_lower_condition_95_observed_boundary, ...
                                  ushift_observed_mean_boundary < ushift_lower_condition_95_observed_boundary, ...
                                  dshift_observed_mean_boundary < dshift_lower_condition_95_observed_boundary];
    upper_condition_95_observed_success = [lshift_observed_mean_boundary >= lshift_upper_condition_95_observed_boundary, ... 
                                  rshift_observed_mean_boundary >= rshift_upper_condition_95_observed_boundary, ...
                                  ushift_observed_mean_boundary >= ushift_upper_condition_95_observed_boundary, ...
                                  dshift_observed_mean_boundary >= dshift_upper_condition_95_observed_boundary];
                              
    % Sampling the cohen d process along the true boundary u = c
    lshift_cohen_d_boundary = abs(lshift_w1.*cohen_d(lshift) + lshift_w2.*cohen_d(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_cohen_d_boundary = abs(rshift_w1.*cohen_d(rshift) + rshift_w2.*cohen_d(rshift(:,[2:dim(2) 1])));
    ushift_cohen_d_boundary = abs(ushift_w1.*cohen_d(ushift) + ushift_w2.*cohen_d(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_cohen_d_boundary = abs(dshift_w1.*cohen_d(dshift) + dshift_w2.*cohen_d(dshift([2:dim(1) 1],:)));
    
    % Calculating the subset condition when residuals in multiplier
    % bootstrap are taken along the cohen d boundary
    lower_condition_80_cohen_d = thr - supGa_cohen_d_80*tau*cohen_d_std;
    upper_condition_80_cohen_d = thr + supGa_cohen_d_80*tau*cohen_d_std;
    lower_condition_90_cohen_d = thr - supGa_cohen_d_90*tau*cohen_d_std;
    upper_condition_90_cohen_d = thr + supGa_cohen_d_90*tau*cohen_d_std;
    lower_condition_95_cohen_d = thr - supGa_cohen_d_95*tau*cohen_d_std;
    upper_condition_95_cohen_d = thr + supGa_cohen_d_95*tau*cohen_d_std;
    
    lshift_lower_condition_80_cohen_d_boundary = abs(lshift_w1.*lower_condition_80_cohen_d(lshift) + lshift_w2.*lower_condition_80_cohen_d(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_lower_condition_80_cohen_d_boundary = abs(rshift_w1.*lower_condition_80_cohen_d(rshift) + rshift_w2.*lower_condition_80_cohen_d(rshift(:,[2:dim(2) 1])));
    ushift_lower_condition_80_cohen_d_boundary = abs(ushift_w1.*lower_condition_80_cohen_d(ushift) + ushift_w2.*lower_condition_80_cohen_d(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_lower_condition_80_cohen_d_boundary = abs(dshift_w1.*lower_condition_80_cohen_d(dshift) + dshift_w2.*lower_condition_80_cohen_d(dshift([2:dim(1) 1],:)));
    
    lshift_upper_condition_80_cohen_d_boundary = abs(lshift_w1.*upper_condition_80_cohen_d(lshift) + lshift_w2.*upper_condition_80_cohen_d(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_upper_condition_80_cohen_d_boundary = abs(rshift_w1.*upper_condition_80_cohen_d(rshift) + rshift_w2.*upper_condition_80_cohen_d(rshift(:,[2:dim(2) 1])));
    ushift_upper_condition_80_cohen_d_boundary = abs(ushift_w1.*upper_condition_80_cohen_d(ushift) + ushift_w2.*upper_condition_80_cohen_d(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_upper_condition_80_cohen_d_boundary = abs(dshift_w1.*upper_condition_80_cohen_d(dshift) + dshift_w2.*upper_condition_80_cohen_d(dshift([2:dim(1) 1],:)));    
    
    lshift_lower_condition_90_cohen_d_boundary = abs(lshift_w1.*lower_condition_90_cohen_d(lshift) + lshift_w2.*lower_condition_90_cohen_d(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_lower_condition_90_cohen_d_boundary = abs(rshift_w1.*lower_condition_90_cohen_d(rshift) + rshift_w2.*lower_condition_90_cohen_d(rshift(:,[2:dim(2) 1])));
    ushift_lower_condition_90_cohen_d_boundary = abs(ushift_w1.*lower_condition_90_cohen_d(ushift) + ushift_w2.*lower_condition_90_cohen_d(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_lower_condition_90_cohen_d_boundary = abs(dshift_w1.*lower_condition_90_cohen_d(dshift) + dshift_w2.*lower_condition_90_cohen_d(dshift([2:dim(1) 1],:)));
    
    lshift_upper_condition_90_cohen_d_boundary = abs(lshift_w1.*upper_condition_90_cohen_d(lshift) + lshift_w2.*upper_condition_90_cohen_d(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_upper_condition_90_cohen_d_boundary = abs(rshift_w1.*upper_condition_90_cohen_d(rshift) + rshift_w2.*upper_condition_90_cohen_d(rshift(:,[2:dim(2) 1])));
    ushift_upper_condition_90_cohen_d_boundary = abs(ushift_w1.*upper_condition_90_cohen_d(ushift) + ushift_w2.*upper_condition_90_cohen_d(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_upper_condition_90_cohen_d_boundary = abs(dshift_w1.*upper_condition_90_cohen_d(dshift) + dshift_w2.*upper_condition_90_cohen_d(dshift([2:dim(1) 1],:)));
    
    lshift_lower_condition_95_cohen_d_boundary = abs(lshift_w1.*lower_condition_95_cohen_d(lshift) + lshift_w2.*lower_condition_95_cohen_d(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_lower_condition_95_cohen_d_boundary = abs(rshift_w1.*lower_condition_95_cohen_d(rshift) + rshift_w2.*lower_condition_95_cohen_d(rshift(:,[2:dim(2) 1])));
    ushift_lower_condition_95_cohen_d_boundary = abs(ushift_w1.*lower_condition_95_cohen_d(ushift) + ushift_w2.*lower_condition_95_cohen_d(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_lower_condition_95_cohen_d_boundary = abs(dshift_w1.*lower_condition_95_cohen_d(dshift) + dshift_w2.*lower_condition_95_cohen_d(dshift([2:dim(1) 1],:)));
    
    lshift_upper_condition_95_cohen_d_boundary = abs(lshift_w1.*upper_condition_95_cohen_d(lshift) + lshift_w2.*upper_condition_95_cohen_d(lshift(:,[dim(2) 1:dim(2)-1])));
    rshift_upper_condition_95_cohen_d_boundary = abs(rshift_w1.*upper_condition_95_cohen_d(rshift) + rshift_w2.*upper_condition_95_cohen_d(rshift(:,[2:dim(2) 1])));
    ushift_upper_condition_95_cohen_d_boundary = abs(ushift_w1.*upper_condition_95_cohen_d(ushift) + ushift_w2.*upper_condition_95_cohen_d(ushift([dim(1) 1:dim(1)-1],:)));
    dshift_upper_condition_95_cohen_d_boundary = abs(dshift_w1.*upper_condition_95_cohen_d(dshift) + dshift_w2.*upper_condition_95_cohen_d(dshift([2:dim(1) 1],:)));
    
    lower_condition_80_cohen_d_success = [lshift_cohen_d_boundary < lshift_lower_condition_80_cohen_d_boundary, ... 
                                  rshift_cohen_d_boundary < rshift_lower_condition_80_cohen_d_boundary, ...
                                  ushift_cohen_d_boundary < ushift_lower_condition_80_cohen_d_boundary, ...
                                  dshift_cohen_d_boundary < dshift_lower_condition_80_cohen_d_boundary];
    upper_condition_80_cohen_d_success = [lshift_cohen_d_boundary >= lshift_upper_condition_80_cohen_d_boundary, ... 
                                  rshift_cohen_d_boundary >= rshift_upper_condition_80_cohen_d_boundary, ...
                                  ushift_cohen_d_boundary >= ushift_upper_condition_80_cohen_d_boundary, ...
                                  dshift_cohen_d_boundary >= dshift_upper_condition_80_cohen_d_boundary];
                              
    lower_condition_90_cohen_d_success = [lshift_cohen_d_boundary < lshift_lower_condition_90_cohen_d_boundary, ... 
                                  rshift_cohen_d_boundary < rshift_lower_condition_90_cohen_d_boundary, ...
                                  ushift_cohen_d_boundary < ushift_lower_condition_90_cohen_d_boundary, ...
                                  dshift_cohen_d_boundary < dshift_lower_condition_90_cohen_d_boundary];
    upper_condition_90_cohen_d_success = [lshift_cohen_d_boundary >= lshift_upper_condition_90_cohen_d_boundary, ... 
                                  rshift_cohen_d_boundary >= rshift_upper_condition_90_cohen_d_boundary, ...
                                  ushift_cohen_d_boundary >= ushift_upper_condition_90_cohen_d_boundary, ...
                                  dshift_cohen_d_boundary >= dshift_upper_condition_90_cohen_d_boundary];
                              
    lower_condition_95_cohen_d_success = [lshift_cohen_d_boundary < lshift_lower_condition_95_cohen_d_boundary, ... 
                                  rshift_cohen_d_boundary < rshift_lower_condition_95_cohen_d_boundary, ...
                                  ushift_cohen_d_boundary < ushift_lower_condition_95_cohen_d_boundary, ...
                                  dshift_cohen_d_boundary < dshift_lower_condition_95_cohen_d_boundary];
    upper_condition_95_cohen_d_success = [lshift_cohen_d_boundary >= lshift_upper_condition_95_cohen_d_boundary, ... 
                                  rshift_cohen_d_boundary >= rshift_upper_condition_95_cohen_d_boundary, ...
                                  ushift_cohen_d_boundary >= ushift_upper_condition_95_cohen_d_boundary, ...
                                  dshift_cohen_d_boundary >= dshift_upper_condition_95_cohen_d_boundary];
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by only comparing
    % binarized sets for residuals on the observed boundary in mult. bootstrap
    if sum(upper_subset_mid_observed_80(:))+sum(mid_subset_lower_observed_80(:))==0
      subset_success_vector_observed_80(t) = 1;
      fprintf('observed nominal 80 success! \n');
    else 
      subset_success_vector_observed_80(t) = 0; 
      fprintf('observed nominal 80 failure! \n');
    end
    
    if sum(upper_subset_mid_observed_90(:))+sum(mid_subset_lower_observed_90(:))==0
      subset_success_vector_observed_90(t) = 1;
      fprintf('observed nominal 90 success! \n');
    else 
      subset_success_vector_observed_90(t) = 0; 
      fprintf('observed nominal 90 failure! \n');
    end

    if sum(upper_subset_mid_observed_95(:))+sum(mid_subset_lower_observed_95(:))==0
      subset_success_vector_observed_95(t) = 1;
      fprintf('observed nominal 95 success! \n');
    else 
      subset_success_vector_observed_95(t) = 0; 
      fprintf('observed nominal 95 failure! \n');
    end
    
    % Testing the subset condition (Ac^- < Ac < Ac^+) by only comparing
    % binarized sets for residuals on the cohen d boundary in mult. bootstrap
    if sum(upper_subset_mid_cohen_d_80(:))+sum(mid_subset_lower_cohen_d_80(:))==0
      subset_success_vector_cohen_d_80(t) = 1;
      fprintf('cohen d nominal 80 success! \n');
    else 
      subset_success_vector_cohen_d_80(t) = 0; 
      fprintf('cohen d nominal 80 failure! \n');
    end
    
    if sum(upper_subset_mid_cohen_d_90(:))+sum(mid_subset_lower_cohen_d_90(:))==0
      subset_success_vector_cohen_d_90(t) = 1;
      fprintf('cohen d nominal 90 success! \n');
    else 
      subset_success_vector_cohen_d_90(t) = 0; 
      fprintf('cohen d nominal 90 failure! \n');
    end

    if sum(upper_subset_mid_cohen_d_95(:))+sum(mid_subset_lower_cohen_d_95(:))==0
      subset_success_vector_cohen_d_95(t) = 1;
      fprintf('cohen d nominal 95 success! \n');
    else 
      subset_success_vector_cohen_d_95(t) = 0; 
      fprintf('cohen d nominal 95 failure! \n');
    end
        
    % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the observed boundary
    if sum(upper_subset_mid_observed_80(:))+sum(mid_subset_lower_observed_80(:)+sum(upper_condition_80_observed_success)+sum(lower_condition_80_observed_success))==0
      subset_success_vector_observed_80_alternate(t) = 1;
      fprintf('observed nominal 80 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_80_alternate(t) = 0; 
      fprintf('observed nominal 80 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_observed_90(:))+sum(mid_subset_lower_observed_90(:)+sum(upper_condition_90_observed_success)+sum(lower_condition_90_observed_success))==0
      subset_success_vector_observed_90_alternate(t) = 1; 
      fprintf('observed nominal 90 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_90_alternate(t) = 0; 
      fprintf('observed nominal 90 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_observed_95(:))+sum(mid_subset_lower_observed_95(:)+sum(upper_condition_95_observed_success)+sum(lower_condition_95_observed_success))==0
      subset_success_vector_observed_95_alternate(t) = 1; 
      fprintf('observed nominal 95 alternate true boundary success! \n');
    else 
      subset_success_vector_observed_95_alternate(t) = 0; 
      fprintf('observed nominal 95 alternate true boundary failure! \n');
    end 
    
        % Testing the subset condition (Ac^- < Ac < Ac^+) by comparing
    % binarized sets as well as the linear interpolated boundary method for
    % residuals taken along the cohen d boundary
    if sum(upper_subset_mid_cohen_d_80(:))+sum(mid_subset_lower_cohen_d_80(:)+sum(upper_condition_80_cohen_d_success)+sum(lower_condition_80_cohen_d_success))==0
      subset_success_vector_cohen_d_80_alternate(t) = 1;
      fprintf('cohen d nominal 80 alternate true boundary success! \n');
    else 
      subset_success_vector_cohen_d_80_alternate(t) = 0; 
      fprintf('cohen d nominal 80 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_cohen_d_90(:))+sum(mid_subset_lower_cohen_d_90(:)+sum(upper_condition_90_cohen_d_success)+sum(lower_condition_90_cohen_d_success))==0
      subset_success_vector_cohen_d_90_alternate(t) = 1; 
      fprintf('cohen d nominal 90 alternate true boundary success! \n');
    else 
      subset_success_vector_cohen_d_90_alternate(t) = 0; 
      fprintf('cohen d nominal 90 alternate true boundary failure! \n');
    end 

    if sum(upper_subset_mid_cohen_d_95(:))+sum(mid_subset_lower_cohen_d_95(:)+sum(upper_condition_95_cohen_d_success)+sum(lower_condition_95_cohen_d_success))==0
      subset_success_vector_cohen_d_95_alternate(t) = 1; 
      fprintf('cohen d nominal 95 alternate true boundary success! \n');
    else 
      subset_success_vector_cohen_d_95_alternate(t) = 0; 
      fprintf('cohen d nominal 95 alternate true boundary failure! \n');
    end                              
end

percentage_success_vector_observed_80                    = mean(subset_success_vector_observed_80, 1);
percentage_success_vector_observed_90                    = mean(subset_success_vector_observed_90, 1);
percentage_success_vector_observed_95                    = mean(subset_success_vector_observed_95, 1);

percentage_success_vector_cohen_d_80                    = mean(subset_success_vector_cohen_d_80, 1);
percentage_success_vector_cohen_d_90                    = mean(subset_success_vector_cohen_d_90, 1);
percentage_success_vector_cohen_d_95                    = mean(subset_success_vector_cohen_d_95, 1);

percentage_success_vector_observed_80_alternate          = mean(subset_success_vector_observed_80_alternate, 1);
percentage_success_vector_observed_90_alternate          = mean(subset_success_vector_observed_90_alternate, 1);
percentage_success_vector_observed_95_alternate          = mean(subset_success_vector_observed_95_alternate, 1);

percentage_success_vector_cohen_d_80_alternate          = mean(subset_success_vector_cohen_d_80_alternate, 1);
percentage_success_vector_cohen_d_90_alternate          = mean(subset_success_vector_cohen_d_90_alternate, 1);
percentage_success_vector_cohen_d_95_alternate          = mean(subset_success_vector_cohen_d_95_alternate, 1);

eval(['save ' SvNm ' nSubj nRlz dim smo mag rimFWHM thr nBoot '... 
      'threshold_cohen_d_80_store threshold_cohen_d_90_store threshold_cohen_d_95_store threshold_observed_80_store threshold_observed_90_store threshold_observed_95_store '...
      'lower_contour_cohen_d_80_store lower_contour_cohen_d_90_store lower_contour_cohen_d_95_store lower_contour_observed_80_store lower_contour_observed_90_store lower_contour_observed_95_store '...
      'upper_contour_cohen_d_80_store upper_contour_cohen_d_90_store upper_contour_cohen_d_95_store upper_contour_observed_80_store upper_contour_observed_90_store upper_contour_observed_95_store '...
      'upper_subset_mid_cohen_d_80_store upper_subset_mid_cohen_d_90_store upper_subset_mid_cohen_d_95_store upper_subset_mid_observed_80_store upper_subset_mid_observed_90_store upper_subset_mid_observed_95_store '...
      'mid_subset_lower_cohen_d_80_store mid_subset_lower_cohen_d_90_store mid_subset_lower_cohen_d_95_store mid_subset_lower_observed_80_store mid_subset_lower_observed_90_store mid_subset_lower_observed_95_store '...
      'subset_success_vector_cohen_d_80 subset_success_vector_cohen_d_90 subset_success_vector_cohen_d_95 subset_success_vector_observed_80 subset_success_vector_observed_90 subset_success_vector_observed_95 subset_success_vector_cohen_d_80_alternate subset_success_vector_cohen_d_90_alternate subset_success_vector_cohen_d_95_alternate subset_success_vector_observed_80_alternate subset_success_vector_observed_90_alternate subset_success_vector_observed_95_alternate '...
      'percentage_success_vector_cohen_d_80 percentage_success_vector_cohen_d_90 percentage_success_vector_cohen_d_95 percentage_success_vector_observed_80 percentage_success_vector_observed_90 percentage_success_vector_observed_95 percentage_success_vector_cohen_d_80_alternate percentage_success_vector_cohen_d_90_alternate percentage_success_vector_cohen_d_95_alternate percentage_success_vector_observed_80_alternate percentage_success_vector_observed_90_alternate percentage_success_vector_observed_95_alternate '...
      'supG_cohen_d_store supG_observed_store '...
      'middle_contour_volume '...
      'lower_contour_cohen_d_80_volume_prct_store lower_contour_cohen_d_90_volume_prct_store lower_contour_cohen_d_95_volume_prct_store lower_contour_observed_80_volume_prct_store lower_contour_observed_90_volume_prct_store lower_contour_observed_95_volume_prct_store '...
      'upper_contour_cohen_d_80_volume_prct_store upper_contour_cohen_d_90_volume_prct_store upper_contour_cohen_d_95_volume_prct_store upper_contour_observed_80_volume_prct_store upper_contour_observed_90_volume_prct_store upper_contour_observed_95_volume_prct_store'])