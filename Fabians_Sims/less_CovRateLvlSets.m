function [coveringRate] = CovRateLvlSets( truef, hatf, thresh, c )

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

truefm = reshape( repmat(reshape( truef, [prod(sf) 1] ), 1, M ) , [sf M] ) ;

coveringRate = 1-sum(any(reshape(...
                ( truefm<c & (hatf >= squeeze(thresh(index{:},2))) ) |...
                ( truefm>=c & (hatf <= squeeze(thresh(index{:},1))))...
               ,[prod(sf) M]), 1 )) / M ;

