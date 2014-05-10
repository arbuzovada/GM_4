function [E, D, M, S] = gibbsIsing(H, J, betaAll, num_iter, connect_type) 
% This function implements Gibbs scheme
%
% INPUT:
%    H: vS-by-vH matrix of double, magnetic field
%    J: 1 or -1, parameter
%    betaAll: 1-by-beta0 array of double, betas
%    num_iter: number of iterations
%    connect_type: 4 or 6, type of neighborhood system
%
% OUTPUT:
%    E: 1-by-beta0 array of double, (math. exp. E(X)) / N
%    D: 1-by-beta0 array of double, (var. E(X)) ^ 0.5 / N
%    M: 1-by-beta0 array of double, (math. exp. mu(X) ^ 2) ^ 0.5 / N
%    S: vS-by-vH-by-beta0 array of ±1, examples of X

    STEP = 20;
    
    [vS, hS] = size(H);
    N = vS * hS;
    beta0 = length(betaAll);
    
    % initialization
    [net, edges] = get_neighbors(vS, hS, beta0, connect_type);
    % CHANGE TO RANDI! ALL ONES FOR SMOOTHNESS
%     S = ones([vS, hS, beta0]);
    S = repmat(2 * randi(2, [vS, hS]) - 3, [1, 1, beta0]);
    nsamples = fix(num_iter / STEP) - fix(num_iter / 3 / STEP);
    E_samples = zeros(nsamples, beta0);
    mu_samples = zeros(nsamples, beta0);
    sample = 0;
    t0 = num_iter - fix(2 * num_iter / 3);
    
    for t = 1 : num_iter
        for i = 1 : N
            % get new distribution
            % p_i = p(x_i = 1 | X_{\i})
            p_i = 1 ./ (1 + exp(-2 * betaAll .* (J * sum(S(net{i}), 1) + H(i))));
            S([0 : beta0 - 1] * N + i) = 2 * (rand(1, beta0) < p_i) - 1;
        end
        if t > t0 && mod(t, STEP) == 0
            sample = sample + 1;
            E_samples(sample, :) = get_E(S, J, H, edges);
            mu_samples(sample, :) = get_mu(S);
        end
    end
    assert(sample == nsamples);
    
    % get E, D, M
    E = mean(E_samples, 1) / N;
    D = (mean(E_samples .^ 2, 1) - mean(E_samples, 1) .^ 2) .^ 0.5 / N;
    M = mean(mu_samples .^ 2, 1) .^ 0.5;
end

function res = get_mu(S)
% This function calculates mu = sum(x_i) / N
    res = squeeze(sum(sum(S)) / (size(S, 1) * size(S, 2)));
end

function res = get_E(S, J, H, edges)
% This function calculates E(X)
    [vS, hS, beta0] = size(S);
    sum_H = squeeze(sum(sum(repmat(H, [1, 1, beta0]) .* S)))';
    S = reshape(S, [vS * hS, beta0]);
    sum_J = sum(S(edges(:, 1), :) .* S(edges(:, 2), :));
    res = -(J * sum_J + sum_H);
end