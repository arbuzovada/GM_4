function [E, D, M, S] = gibbsIsing4(H, J, betaAll, num_iter) 
% This function implements Gibbs scheme
%
% INPUT:
%    H: vS-by-vH matrix of double, magnetic field
%    J: 1 or -1, parameter
%    betaAll: 1-by-beta0 array of double, betas
%    num_iter: number of iterations
%
% OUTPUT:
%    E: 1-by-beta0 array of double, (math. exp. E(X)) / N
%    D: 1-by-beta0 array of double, (var. E(X)) ^ 0.5 / N
%    M: 1-by-beta0 array of double, (math. exp. mu(X) ^ 2) ^ 0.5 / N
%    S: vS-by-vH-by-beta0 array of ±1, examples of X

    [vS, vH] = size(H);
    N = vS * vH;
    beta0 = length(betaAll);
    
    % initialization
    E = zeros(1, beta0);
    D = zeros(1, beta0);
    M = zeros(1, beta0);
    S = repmat(2 * randi(2, [vS, vH]) - 3, [1, 1, beta0]);
    E_samples = zeros(fix(2 * num_iter / 3), beta0);
    mu_samples = zeros(fix(2 * num_iter / 3), beta0);
    t0 = num_iter - fix(2 * num_iter / 3);
    
    for t = 1 : num_iter
        for i = 1 : N
            % get new distribution
            % p_i = p(x_i = 1 | X_{\i})
            p_i = 1 ./ (1 + exp(-2 * betaAll .* (J * ngbrs(i, S) + H(i))));
            S([0 : beta0 - 1] * N + i) = 2 * (rand(1, beta0) < p_i) - 1;
        end
        if t > t0
            E_samples(t - t0, :) = get_E(S, J, H);
            mu_samples(t - t0, :) = get_mu(S);
        end
    end
    
    % get E, D, M
    E = mean(E_samples, 1) / N;
    D = mean((E_samples - repmat(E, size(E_samples, 1), 1)) .^ 2) .^ 0.5 / N;
    M = mean(mu_samples .^ 2, 1) .^ 0.5;
end

function res = ngbrs(i, S)
% This function returns sum of values of x_i neighbors for each beta
    [vS, vH, beta0] = size(S);
    ngbr = [i + 1, i - 1, i + vS, i - vS];
    mask = (ngbr > 0) & (ngbr <= vS * vH);
    mask(1) = mask(1) & (mod(i, vS) ~= 0);
    mask(2) = mask(2) & (mod(i - 1, vS) ~= 0);
    n = sum(mask);
    inds = repmat([0 : beta0 - 1]' * vS * vH, 1, n) + ...
        repmat(ngbr(mask), beta0, 1);
    res = sum(reshape(S(inds), [beta0, n]), 2)';
end

function res = get_mu(S)
% This function calculates mu = sum(x_i) / N
    res = squeeze(sum(sum(S)) / (size(S, 1) * size(S, 2)));
end

function res = get_E(S, J, H)
% This function calculates E(X)
    [vS, vH, beta0] = size(S);
    N = vS * vH;
    res = zeros(1, beta0);
    for i = 1 : N
        cur = 0;
        if (mod(i, vS) ~= 0) % down
            cur = cur + S([0 : beta0 - 1] * N + i + 1);
        end
        if (fix(i / vS) < vH - 1) % right
            cur = cur + S([0 : beta0 - 1] * N + i + vS);
        end
        res = res + S([0 : beta0 - 1] * N + i) .* (J * cur + H(i));
    end
    res = -res;
end