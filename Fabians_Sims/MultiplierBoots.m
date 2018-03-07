function [quantiles] = MultiplierBoots( R, alpha, Mboot, mask, center, normalize )

% Bootstraps the alpha-quantile of the maximum of a Gaussian field
% Input:
% R:    random field over a domain in R^D, it is an (D+1)-dimensional array,
%       where the last dimension enumerates the samples
% alpha:     vector for which the quantile needs to be estimated
% Mboot:     amount of bootstrap replicates (default=1e3)
% mask:      array containing the values which are inside the domain (default=ones(size(R)(1:2)))
% center:    option to center the field using the sample mean (default=1)
% normalize: option to normalize the field by sample variance (default=1)
%Output:
% quantile is the bootstrapped quantile of the maximum distribution of the 
% input processes

% Check number of inputs.
if nargin > 6
    error('MultiplierBoots requires at most 4 optional inputs');
end

% Fill in unset optional values.
switch nargin
    case 2
        Mboot     = 5e3;
        mask      = ones(size(R,1),size(R,2));
        center    = 1;
        normalize = 1;
    case 3
        mask      = ones(size(R,1),size(R,2));
        center    = 1;
        normalize = 1;
    case 4
        center    = 1;
        normalize = 1;
    case 5
        normalize = 1;
end

%%%%%% Get parameters for simulation
dimR = size(R) ;
N    = dimR(end) ;

% Combine values of R in the mask to an matrix for faster computation
R       = reshape(R, prod(dimR(1:end-1)), N);
vecmask = repmat( reshape(logical(mask), prod(dimR(1:end-1)), 1), [1 N] );

R = R(vecmask) ;
R = reshape(R, [length(R)/N N] );
      
% Center/normalize, if required
if center
    R = R - repmat( mean(R,2), [1 N] );
end
if normalize
   R = R ./ repmat( transpose(sqrt(var( transpose(R)))), [1 N] );
end
    

%%%%% compute bootstrap replicates
% compute multipliers
multiplier = normrnd( 0, 1, [N,Mboot] );

% compute bootstrapped means
bootMeans      = R * multiplier / N;
bootSecMoments = R.^2 * multiplier.^2 / N;
% we put an abs here to make sure that no NaNs are produced due to machine precision error.
bootSigma = sqrt( abs(bootSecMoments - bootMeans.^2) / (N-1) * N );

%compute bootstrapped values of maximum
bootMax = max(abs( sqrt(N)*bootMeans./bootSigma ));

%compute quantiles from the bootstrapp distribution of maximum
quantiles = quantile( bootMax, alpha );