function [coveringRate] = AB_CovRateLvlSets( truef, hatf, thresh, c )

% Compute Euler characteristic of excursion set
% Input:
% truef:  underlying true function on an N-dimensional array, it is saved as
%         an N-dimensional array, where the values correspond to the heights
% hatf:   array of size [size(truef), M] containing M estimates of truef
% thresh: array of size [size(truef), M, 2] corresponding to the thresholds for the
%         hatf, [ , , 1] are lower bounds, [ , , 2] are upper bounds
% c:      targeted levelset
%Output:
% coveringRate is computed from all available hatf and computes the
% frequency that the true c-levelset of truef is completly contained in the
% thresholded hatfs.

sf    = size( truef );
N     = ndims( truef );
M     = size( hatf, N+1 );
index = repmat( {':'}, 1, N+1 );
middle_contour = truef >= c;
subset_success_vector = zeros(M,1);

for i=1:M
    lower_contour = squeeze(hatf(:,:,i)) >= squeeze(thresh(:,:,i,1));
    upper_contour = squeeze(hatf(:,:,i)) >= squeeze(thresh(:,:,i,2));
    
    mid_on_upper              = upper_contour.*middle_contour;
    lower_on_mid              = middle_contour.*lower_contour;
    upper_subset_mid          = upper_contour - mid_on_upper;
    mid_subset_lower          = middle_contour - lower_on_mid;
    
    if sum(upper_subset_mid(:))+sum(mid_subset_lower(:))==0
      subset_success_vector(i) = 1; 
    else 
      subset_success_vector(i) = 0;
    end
end

coveringRate = mean(subset_success_vector, 1);
