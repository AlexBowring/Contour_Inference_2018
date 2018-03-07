%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%        Covering rate of CoPE sets
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

% Set working path
cd '/storage/maullz/Contour_Inference_2018/Fabians_Sims'

% Field parameters
NOISE_TYPE = 'isotropic';
% NOISE_TYPE = 'scale-space';
L0 = 1;     % EC of domain
D = 2;      % Dimension of domain (including scale)
L = 51;    % size of domain
N = 500;     % replicates (sample size)
nsim = 1000;  % number of simulations
nu = 5;     % parameter for isotropic
gamma = 4:.2:40;  % bandwiths for scale-space
b = 5;      % size of filter domain in std
sigma = 1;

% Signal
SIGNAL_TYPE = 'signal'%'SNR';
SIGNAL_SHAPE = 'linear';
aa = 4;      % Slope of signal


c = 2;      % Target SNR level


%% Generate fields
switch(NOISE_TYPE),
    case 'isotropic', params = nu;
    case 'scale-space', params = gamma;
end

% Noise
tic
[eps, LKC] = generateField(NOISE_TYPE, N, nsim, D, L, params, b);
toc
% Save error processes for later use (files might get large (N=200, nsim=1000 => 4gB!))
save('N200Nsim1000_isotropic_nu5.mat', '-v7.3')
%% Compute Signal and Observations
t = (0.5:L)/L;
[xx, yy] = meshgrid(t, t);
switch(SIGNAL_SHAPE),
    case 'linear', mu = aa*xx;
    case 'quadratic', mu = aa*(1 - xx.^2 - yy.^2);
end
y = repmat(mu, [1 1 N nsim]) + eps;

% Observed
switch(SIGNAL_TYPE),
    case 'signal', delta = mu; deltahat = mean(y,3);
    case 'SNR', delta = mu/sigma; deltahat = mean(y, 3) ./ std(y, 0, 3);
end
deltahat = permute(deltahat, [1 2 4 3]);
clear xx yy t

%% Simulation of covering rate
% Evaluate maximum at the correct boundary
mask = (delta >= c) & (delta <= c);

% Initialize vector for threshold a
a = zeros(1,nsim);

%%%%% Simulation of estimated quantiles
tic
for k = 1:nsim
    % Obtain quantile estimate
    a(k) = MultiplierBoots( y(:,:,:,k), 0.95, 5000, mask, 1, 1 );
end
toc
%%
% estimate empirical variance of the fields
hatsigma = permute(sqrt(var(y, 1, 3)), [1 2 4 3]);
% compute upper and lower threshold for CoPE sets
thresh = ones([size(delta) nsim 2]);
thresh(:,:,:,1) = c - repmat(shiftdim(a, -1), [51 51 1]).*hatsigma / sqrt(N);
thresh(:,:,:,2) = c + repmat(shiftdim(a, -1), [51 51 1]).*hatsigma / sqrt(N);
% Compute the covering rate
[covRate]            = CovRateLvlSets( delta, deltahat, thresh, c )
[lessthan_covRate]   = less_CovRateLvlSets( delta, deltahat, thresh, c )
[ABcovRate]          = AB_CovRateLvlSets( delta, deltahat, thresh, c )

%% Figures
ll=8 % pick 11th simulation for plots
figure(1)
imagesc(y(:,:,1)), colorbar
title('Observed field')

figure(2)
imagesc(delta), colorbar
title('Delta')

figure(3), hold on
imagesc(deltahat(:,:,ll)), colorbar
contour(delta, [1 1]*c, 'r', 'Linewidth', 2)
contour(deltahat(:,:,ll), [1 1]*c, 'k', 'Linewidth', 2)
contour(deltahat(:,:,ll)-thresh(:,:,ll,1), [1 1]*0, 'k--', 'Linewidth', 2)
contour(deltahat(:,:,ll)-thresh(:,:,ll,2), [1 1]*0, 'k--', 'Linewidth', 2)
title('Excursion sets'), hold off

%%
figure(5)
imagesc(delta<=c & (deltahat(:,:,1) >= thresh(:,:,ll,2)))
figure(6)
imagesc(delta>=c & (deltahat(:,:,1) <= thresh(:,:,ll,1)))

figure(7)
imagesc( delta<=c & (deltahat(:,:,1) >= thresh(:,:,ll,2)) | delta>=c & (deltahat(:,:,1) <= thresh(:,:,ll,1)) )
title('Locations, where covering is invalid')
