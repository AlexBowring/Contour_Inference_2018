function [f, LKC] = generateField(TYPE, n, nsim, N, L, params, b, L0)

% Generate random fields
% Input:
%   TYPE:   Type of field ('isotropic' or 'scale-space')
%   n:      Sample size
%   nsim:   Number of simulations
%   N:      Dimension of domain (including scale)
%   L:      Size of field domain (default = 200)
%   params: nu (scalar) for isotropic, gamma (vector) for scale-space
%   b:      half of size of filter domain in std (default = 5)
%   L0:     Euler characteristic of domain (for scale-space) (default = 1)
%
% Output:
%   f:      Array of size Lx()xnxnsim
%   LKC:    Theoretical LKCs

% Check input
if ~strcmp(TYPE, 'isotropic') & ~strcmp(TYPE, 'scale-space'),
    error([TYPE, ' not implemented']),
end
if ~exist('N'), N = 1; end
if strcmp(TYPE, 'isotropic') & (N ~= 1 & N ~=2 & N~=3), error('N must be 1 or 2'); end
if strcmp(TYPE, 'scale-space') & (N ~= 2), error('N must be 2'); end
if ~exist('L'), L = 50; end
if ~exist('params') & strcmp(TYPE, 'isotropic'), params = 5; end
if ~exist('params') & strcmp(TYPE, 'scale-space'), params = 4:.2:40; end
if ~exist('b'), b = 5; end
if ~exist('L0'), L0 = 1; end

% Computations
switch TYPE,
case 'isotropic',
    % LKCs
    nu = params;
    sigma_nu = ( 2^N * nu^N * pi^(N/2) )^(-1/2); % SD of the smoothed field
    alpha = 1/(4*nu^2);
    switch N,
    case 1,
        LKC = (L-1) * sqrt(2*alpha) ;
    case 2,
        LKC = [2*(L-1); (L-1)^2] .* [sqrt(2*alpha); 2*alpha] ;
    case 3,
        LKC = [3*(L-1)/2; 3*(L-1)^2; (L-1)^3] .* [sqrt(2*alpha); 2*alpha; (2*alpha)^(3/2)] ;
    end
    
    % Generate fields
    switch N,
        case 1,
            x = -b*nu : b*nu;
            h = (2*pi*nu^2)^(-1/2) * exp(-x.^2 / (2*nu^2));
            % h = fspecial('gaussian', [2*b*nu+1 1], nu); % Equivalent
        case 2,
            h = fspecial('gaussian', (2*b*nu+1)*[1 1], nu);
        case 3,
            h1 = fspecial('gaussian', [2*b*nu+1 1], nu);
            h2 = fspecial('gaussian', (2*b*nu+1)*[1 1], nu);
            len = length(h1);
            h = repmat(h2, [1 1 len]) .* repmat(shiftdim(h1, -2), [len len 1]);
    end
    w = randn([repmat(L + 2*b*nu, 1, N), n, nsim]);
    f = convn(w, h, 'valid');
    f = f / sigma_nu; % standardized field
    
    % f = f / sqrt(sum(h(:).^2)); % standardized field
    % For N=2, sum(h(:).^2) is about 1/(4*pi*nu^2)

case 'scale-space',
    % LKCs (Siegmund 1995)
    NN = N - 1; % Dimension of domain (not including scale)
    kappa = 0.5;
    lambda = 0.5;
    
    % kappa and lambda can be computed as below.
    % du = 0.02;
    % u = -4:du:4;
    % k = pi^(-1/4)*exp(-u.^ 2 / 2);
    % dk = pi^(-1/4)*(-u).*exp(-u.^ 2 / 2);
    % kappa = sum((u.*dk + NN*k/2).^2 .* du);
    % lambda = sum((dk.*dk) .* du);

    gamma = params;
    L2 = (L-1) * (gamma(1)^(-1) - gamma(end)^(-1)) * sqrt(lambda*kappa);
    L1 = (L-1)/2 * (gamma(1)^(-NN) + gamma(end)^(-NN)) * sqrt(lambda) + L0 * log(gamma(end)/gamma(1)) * sqrt(kappa);
    LKC = [L1; L2] ;

    % Generate fields
    whitenoise_large = randn(L + 2*b*gamma(end), n, nsim);
    f = zeros(L, length(gamma), n, nsim);
    for k = 1:length(gamma),
        whitenoise = whitenoise_large((1+round(b*gamma(end)-b*gamma(k))):(L+round(b*gamma(end)+b*gamma(k))), :, :);
        x = (-b*gamma(k):b*gamma(k))';
        w = gamma(k)^(-1/2)*pi^(-1/4)*exp(-x.^ 2 / (2 * gamma(k)^2));
        z_gamma = convn(whitenoise, w, 'valid');
        f(:,k,:,:) = z_gamma;
    end
end



