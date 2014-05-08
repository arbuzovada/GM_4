function [net, edges] = get_neighbors(vS, hS, beta0, connect_type)
% This function find neighbor indices for each vertice
% res: 1-by-(vS * vH) cell array of arrays of indices

    N = vS * hS;
    net = cell(1, N);
    edges = zeros(0, 2);
    for i = 1 : N
        if (mod(i, vS) ~= 0) % down
            net{i} = [net{i}; i + 1];
            edges = [edges; i, i + 1]; %#ok<AGROW>
        end
        if (mod(i - 1, vS) ~= 0) % up
            net{i} = [net{i}; i - 1];
        end
        if (i - vS > 0) % left
            net{i} = [net{i}; i - vS];
        end
        if (i + vS <= N) % right
            net{i} = [net{i}; i + vS];
            edges = [edges; i, i + vS]; %#ok<AGROW>
        end
        if connect_type == 6
            if (mod(i - 1, vS) ~= 0) && (i - vS > 0) % up-left
                net{i} = [net{i}; i - vS - 1];
            end
            if (mod(i, vS) ~= 0) && (i + vS <= N) % down-right
                net{i} = [net{i}; i + vS + 1];
                edges = [edges; i, i + vS + 1]; %#ok<AGROW>
            end
        end
        net{i} = repmat([0 : beta0 - 1] * vS * hS, length(net{i}), 1) + ...
            repmat(net{i}, 1, beta0);
    end
end