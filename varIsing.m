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

    tic
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
    net = get_neighbors(vS, hS, beta0, connect_type);
    % probability of x = -1 and 1 for each temperature
    q = repmat(rand(vS, hS), [1, 1, beta0, 2]);
    % zero represents -1, one represents 1
    integr = cell(1, 6);
    for i = 1 : 6
        integr{i}.sum = sum(2 * de2bi(0 : 2 ^ i - 1) - 1, 2)';
        integr{i}.var = zeros(i, beta0, 2 ^ i);
        integr{i}.var(:, 1, :) = N * beta0 * de2bi(0 : 2 ^ i - 1)';
        integr{i}.var = repmat(integr{i}.var(:, 1, :), [1, beta0, 1]);
    end
    E = zeros(1, beta0);
    D = zeros(1, beta0);
    M = zeros(1, beta0);
    L = zeros(1, beta0);
    toc
    
    for t = 1 : MAX_ITER
        for i = 1 : N
            % get new distribution q_i(x_i)
            m = size(net{i}, 1); % number of neighbors
            cur = zeros(2 ^ m, beta0);
            for j = 1 : 2 ^ m
                cur(j, :) = prod(q(net{i} + integr{m}.var(:, :, j)), 1);
            end
            % x = 1
            q(N * beta0 + [0 : beta0 - 1] * N + i) = 1 ./ ...
                (1 + exp(-2 * betaAll .* ((J * integr{m}.sum + H(m)) * cur)));
            % x = -1
            q([0 : beta0 - 1] * N + i) = 1 - q(N * beta0 + [0 : beta0 - 1] * N + i);
        end
    end
end      