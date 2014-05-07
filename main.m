vS = 20;
hS = 20;
T = 0.1 : 0.1 : 10;
betaAll = 1 ./ T;
H = zeros(vS, hS);%2 * rand(vS, hS) - 1;
J = 1;
num_iter = 100;
opt_params.max_iter = 100;

tic
% [E, D, M, S] = gibbsIsing(H, J, betaAll, num_iter, 4);
[E, D, M, L] = varIsing(H, J, betaAll, opt_params, 4);
toc

% figure()
% plot(T, E)
% hold on
% plot(T, D, 'r')
% plot(T, M, 'g')
% legend('E[X]/N', 'sqrt(D[X])/N', 'sqrt(E[M^2])')
% xlabel('T')