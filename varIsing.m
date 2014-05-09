function [E, D, M, L] = varIsing(H, J, betaAll, opt_params, connect_type)
% This function implements variational inference
%
% INPUT:
%    H: vS-by-vH matrix of double, magnetic field
%    J: 1 or -1, parameter
%    betaAll: 1-by-beta0 array of double, betas
%    opt_params: structure
%        max_iter: number of maximum iterations (default: 300)
%        tol_crit: tolerance in stopping criteria (default: 1e-4)
%        num_start: number of starts (default: 1)
%    connect_type: 4 or 6, type of neighborhood system
%
% OUTPUT:
%    E: 1-by-beta0 array of double, (math. exp. E(X)) / N
%    D: 1-by-beta0 array of double, (var. E(X)) ^ 0.5 / N
%    M: 1-by-beta0 array of double, (math. exp. mu(X) ^ 2) ^ 0.5 / N
%    L: 1-by-beta0 array of double, lower estimates of normalization
%        constant

    MAX_ITER = 300;
    TOL_CRIT = 1e-4;
    NUM_START = 1;
    if nargin > 3
        if isfield(opt_params, 'max_iter')
            MAX_ITER = opt_params.max_iter;
        end
        if isfield(opt_params, 'tol_crit')
            TOL_CRIT = opt_params.tol_crit;
        end
        if isfield(opt_params, 'num_start')
            NUM_START = opt_params.num_start;
        end
    end
    [vS, hS] = size(H);
    N = vS * hS;
    beta0 = length(betaAll);

    % initialization
    [net, edges] = get_neighbors(vS, hS, beta0, connect_type);
    % probability of x = -1 and 1 for each temperature
    q = repmat(rand(N, 1), 1, beta0);
    H = reshape(H, N, 1);
%     q_init = q;
%     load('q_init.mat');
%     q = q_init;
    D = zeros(1, beta0);
    
    for t = 1 : MAX_ITER
        for i = 1 : N 
            % get new distribution q_i(x_i = 1)
            q(i, :) = 1 ./ ...
                (1 + exp(-2 * betaAll .* (J * sum(2 * q(net{i}) - 1, 1) + H(i))));
        end
        % L(q)
        L = betaAll .* (sum(repmat(H, 1, beta0) .* (2 * q - 1)) + ...
            J * sum((2 * q(edges(:, 1), :) - 1) .* (2 * q(edges(:, 2), :) - 1)));
        Llog = log(q) .* q + log(1 - q) .* (1 - q);
        Llog(isnan(Llog)) = 0;
        L = L - sum(Llog);
        if all(L < TOL_CRIT)
            break;
        end
    end
    E = -(sum(repmat(H, 1, beta0) .* (2 * q - 1)) + ...
        J * sum((2 * q(edges(:, 1), :) - 1) .* (2 * q(edges(:, 2), :) - 1))) / N;
    inds = [reshape(repmat([1 : N]', 1, N)', N ^ 2, 1), repmat([1 : N], 1, N)'];
    inds = inds(inds(:, 1) ~= inds(:, 2), :);
    M = (sum((2 * q(inds(:, 1), :) - 1) .* (2 * q(inds(:, 2), :) - 1)) + N) .^ 0.5 / N ^ 2;
end      