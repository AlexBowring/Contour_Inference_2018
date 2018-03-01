function Sim_07(nSubj,SvNm,nRlz)
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
dim     = [51 51]; 
mag     = 3;
smo     = 10;
rimFWHM = 15; 				 
stdblk  = prod(dim([1 2]));
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
Sig = repmat(linspace(1, 3, dim(2)), dim(2), 1);

% Uncomment to look at the Signal
%imagesc(Sig); axis image; colorbar
AC = Sig >= thr;

% The boundary for Sig > 2, note that Sig = 2.02 in the 51st column
true_boundary = zeros(dim);
true_boundary(:,26) = ones(51, 1);
true_boundary = logical(true_boundary);


mean_std_store = zeros(nRlz, 1);
tNoises_store = zeros([dim nSubj]);

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
        tNoises_store(:,:,i) = tNoises;
        tImgs = Sig + tNoises; % Creates the true image of smoothed signal + smoothed noise
        observed_data(:,:,i) = tImgs;
        
      end %========== Loop i (subjects)
      reshape_tNoises = reshape(tNoises_store, [prod(dim) nSubj]);
      mean_std_store(t, 1) = mean(std(reshape_tNoises(:,[1:nSubj])));
end