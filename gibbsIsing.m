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

    [vS, hS] = size(H);
    N = vS * hS;
    beta0 = length(betaAll);
    
    % initialization
    net = get_neighbors(vS, hS, beta0, connect_type);
    % CHANGE TO RANDI! ALL ONES FOR SMOOTHNESS
    S = ones([vS, hS, beta0]); % repmat(2 * randi(2, [vS, hS]) - 3, [1, 1, beta0]);
    E_samples = zeros(fix(2 * num_iter / 3), beta0);
    mu_samples = zeros(fix(2 * num_iter / 3), beta0);
    t0 = num_iter - fix(2 * num_iter / 3);
    
    for t = 1 : num_iter
        for i = 1 : N
            % get new distribution
            % p_i = p(x_i = 1 | X_{\i})
            p_i = 1 ./ (1 + exp(-2 * betaAll .* (J * sum(S(net{i}), 1) + H(i))));
            S([0 : beta0 - 1] * N + i) = 2 * (rand(1, beta0) < p_i) - 1;
        end
        if t > t0
            E_samples(t - t0, :) = get_E(S, J, H, connect_type);
            mu_samples(t - t0, :) = get_mu(S);
        end
    end
    
    % get E, D, M
    E = mean(E_samples, 1) / N;
    D = (mean(E_samples .^ 2, 1) - mean(E_samples, 1) .^ 2) .^ 0.5 / N;
    M = mean(mu_samples .^ 2, 1) .^ 0.5;
end

function res = get_mu(S)
% This function calculates mu = sum(x_i) / N
    res = squeeze(sum(sum(S)) / (size(S, 1) * size(S, 2)));
end

function res = get_E(S, J, H, connect_type)
% This function calculates E(X)
    [~, ~, beta0] = size(S);
    sum_H = squeeze(sum(sum(repmat(H, [1, 1, beta0]) .* S)))';
    if connect_type == 4
        sum_J = squeeze(sum(sum(S(1 : (end - 1), :, :) .* S(2 : end, :, :))) + ...
            sum(sum(S(:, 1 : (end - 1), :) .* S(:, 2 : end, :))))';
    else
        sum_J = squeeze(sum(sum(S(1 : (end - 1), :, :) .* S(2 : end, :, :))) + ...
            sum(sum(S(:, 1 : (end - 1), :) .* S(:, 2 : end, :))) + ...
            sum(sum(S(1 : (end - 1), 1 : (end - 1), :) .* S(2 : end, 2 : end, :))))';
    end
    res = -(J * sum_J + sum_H);
end