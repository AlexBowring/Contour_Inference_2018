function [X,ss] = spm_conv(X,sx,sy)
% Gaussian convolution
% FORMAT [X,ss] = spm_conv(X,sx[,sy]);
% X    - matrix
% sx   - kernel width (FWHM) in pixels
% sy   - optional non-isomorphic smoothing
%
% ss   - sum of squares of kernel
%___________________________________________________________________________
%
% spm_conv is a one or two dimensional convolution of a matrix variable
% in working memory.  It capitalizes on the sparisity structure of the
% problem and the separablity of multidimensional convolution with a Gaussian
% kernel by using one-dimensional convolutions and kernels that are
% restricted to non near-zero values
%
% If image is mean zero and variance one, scaling the smoothed image by
% 1/sqrt(ss) will return it to unit variance.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_conv.m 4195 2011-02-05 18:39:13Z ged $


% assume isomorphic smoothing
%---------------------------------------------------------------------------
if nargin < 3; sy = sx; end
sx      = abs(sx);
sy      = abs(sy);
[lx ly] = size(X);



% FWHM -> sigma
%---------------------------------------------------------------------------
sx    = sx/sqrt(8*log(2)) + eps;
sy    = sy/sqrt(8*log(2)) + eps;

% kernels
%---------------------------------------------------------------------------
Ex    = min([round(6*sx) lx]);
x     = [-Ex:Ex];
kx    = exp(-x.^2/(2*sx^2));
kx    = kx/sum(kx);
Ey    = min([round(6*sy) ly]);
y     = [-Ey:Ey];
ky    = exp(-y.^2/(2*sy^2));
ky    = ky/sum(ky);

[Sx Sy]  = meshgrid(kx,ky);
ss       = sum((Sx(:).*Sy(:)).^2); 

% convolve
%---------------------------------------------------------------------------
if lx > 1;
    for i = 1:ly
        u      = X(:,i);
        u      = [flipud(u(1:Ex)); u; flipud(u([1:Ex] + lx - Ex))];
        U      = sparse(conv(full(u),kx));
        X(:,i) = U([1:lx] + 2*Ex);
    end
end
if ly > 1;
    for i = 1:lx
        u      = X(i,:);
        u      = [fliplr(u(1:Ey)) u fliplr(u([1:Ey] + ly - Ey))];
        U      = sparse(conv(full(u),ky));
        X(i,:) = U([1:ly] + 2*Ey);
    end
end
