function [best_E, best_D, best_M, best_L] = varIsing(H, J, betaAll, ...
    opt_params, connect_type)
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
    H = reshape(H, N, 1);
    best_L = zeros(1, beta0);
    best_q = zeros(N, beta0);
    
    for w = 1 : NUM_START
        % probability of x = -1 and 1 for each temperature
        q = repmat(rand(N, 1), 1, beta0);
        old_L = zeros(1, beta0);
        for t = 1 : MAX_ITER
            for i = 1 : N 
                % get new distribution q_i(x_i = 1)
                q(i, :) = 1 ./ ...
                    (1 + exp(-2 * betaAll .* ...
                    (J * sum(2 * q(net{i}) - 1, 1) + H(i))));
            end
            % L(q)
            L = betaAll .* (J * sum((2 * q(edges(:, 1), :) - 1) .* ...
                (2 * q(edges(:, 2), :) - 1)) + ...
                sum(bsxfun(@times, 2 * q - 1, H)));
            Llog = log(q) .* q + log(1 - q) .* (1 - q);
            Llog(isnan(Llog)) = 0;
            L = L - sum(Llog);
            if all(abs(L - old_L) < TOL_CRIT)
                break;
            end
            old_L = L;
        end
        E = -(sum(bsxfun(@times, 2 * q - 1, H)) + ...
            J * sum((2 * q(edges(:, 1), :) - 1) .* ...
            (2 * q(edges(:, 2), :) - 1))) / N;
        M = (sum(2 * q - 1) .^ 2 - sum((2 * q - 1) .^ 2) + N) .^ 0.5 / N;
        D1 = sum(bsxfun(@times, (2 * q - 1), H)) .^ 2 - ...
            sum(bsxfun(@times, (2 * q - 1), H) .^ 2) + ...
            sum(sum(H .^ 2));
        nedges = size(edges, 1);
        inds = [repmat(edges, N, 1), ...
            reshape(repmat([1 : N]', 1, nedges)', N * nedges, 1)];
        pair_inds = inds(inds(:, 2) == inds(:, 3), :);
        pair_inds = [pair_inds; inds(inds(:, 1) == inds(:, 3), [2 1 3])]; %#ok
        D2 = 2 * J * (sum((2 * q(edges(:, 1), :) - 1) .* ...
            (2 * q(edges(:, 2), :) - 1)) .* ...
            sum(bsxfun(@times, 2 * q - 1, H)) - ...
            sum(bsxfun(@times, ...
            (2 * q(pair_inds(:, 1), :) - 1) .* ...
            (2 * q(pair_inds(:, 2), :) - 1) .* ...
            (2 * q(pair_inds(:, 3), :) - 1), H(pair_inds(:, 3)))) + ...
            sum(bsxfun(@times, ...
            (2 * q(pair_inds(:, 1), :) - 1), H(pair_inds(:, 3)))));
        inds = [reshape(repmat(edges, 1, nedges)', 2, nedges ^ 2)', ...
            repmat(edges', 1, nedges)'];
        % equal edges -> + nedges in D3
        inds = inds(~all((inds(:, 1:2) == inds(:, 3:4))'), :);
        % 1 equal vertice -> pair_inds
        pair_inds = inds(inds(:, 1) == inds(:, 3), :);
        pair_inds = [pair_inds; inds(inds(:, 1) == inds(:, 4), [1 2 4 3])]; %#ok
        pair_inds = [pair_inds; inds(inds(:, 2) == inds(:, 3), [2 1 3 4])]; %#ok
        pair_inds = [pair_inds; inds(inds(:, 2) == inds(:, 4), [2 1 4 3])]; %#ok
        D3 = J ^ 2 * (sum((2 * q(edges(:, 1), :) - 1) .* ... % all possible
            (2 * q(edges(:, 2), :) - 1)) .^ 2 - ...
            sum(((2 * q(edges(:, 1), :) - 1) .* ... % equal edges
            (2 * q(edges(:, 2), :) - 1)) .^ 2) + nedges - ...
            sum((2 * q(pair_inds(:, 1), :) - 1) .^ 2 .* ... % 1 vertice
            (2 * q(pair_inds(:, 2), :) - 1) .* ...
            (2 * q(pair_inds(:, 4), :) - 1)) + ...
            sum((2 * q(pair_inds(:, 2), :) - 1) .* ...
            (2 * q(pair_inds(:, 4), :) - 1)));
        D = (D1 + D2 + D3 - (N * E) .^ 2) .^ 0.5 / N;
        best_inds = L > best_L;
        best_E(best_inds) = E(best_inds);
        best_D(best_inds) = D(best_inds);
        best_M(best_inds) = M(best_inds);
        best_L(best_inds) = L(best_inds);
        best_q(:, best_inds) = q(:, best_inds);
    end
    save q_final.mat best_q vS hS beta0
end