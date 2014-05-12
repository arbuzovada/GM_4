vS = 20;
hS = 20;
T = 0.5 : 0.1 : 10;
betaAll = 1 ./ T;
H = zeros(vS, hS);
J = -1;
num_iter = 10000;
opt_params.max_iter = 500;
opt_params.num_start = 200;

tic
[E4, D4, M4, L4] = varIsing(H, J, betaAll, opt_params, 4);
toc
tic
[E6, D6, M6, L6] = varIsing(H, J, betaAll, opt_params, 6);
toc

figure()
hold on
plot(T, D4, 'b')
plot(T, D6, 'r')
legend('rectangular', 'triangular')
xlabel('$T$', 'interpreter', 'latex')
ylabel('$\frac{1}{N}\sqrt{\mathbf{D}[X]}$', 'interpreter', 'latex');%, 'FontSize', 13
print('varD_m1', '-depsc2', '-r300');
% eps2xxx('varD_1.eps', {'jpeg'}, 'C:\Program Files\gs\gs9.10\bin\gswin64c.exe')

figure()
hold on
plot(T, E4, 'b')
plot(T, E6, 'r')
legend('rectangular', 'triangular')%, 'Location', 'SouthEast')
xlabel('$T$', 'interpreter', 'latex')
ylabel('$\frac{1}{N}\mathbf{E}[X]$', 'interpreter', 'latex');
print('varE_m1', '-depsc2', '-r300');
% eps2xxx('varE_1.eps', {'jpeg'}, 'C:\Program Files\gs\gs9.10\bin\gswin64c.exe')

figure()
hold on
plot(T, M4, 'b')
plot(T, M6, 'r')
legend('rectangular', 'triangular')%, 'Location', 'SouthEast')
xlabel('$T$', 'interpreter', 'latex')
ylabel('$\sqrt{\mathbf{E}[\mu^2(X)]}$', 'interpreter', 'latex');
print('varM_m1', '-depsc2', '-r300');
% eps2xxx('varM_1.eps', {'jpeg'}, 'C:\Program Files\gs\gs9.10\bin\gswin64c.exe')

figure()
hold on
plot(T, L4, 'b')
plot(T, L6, 'r')
legend('rectangular', 'triangular')%, 'Location', 'SouthEast')
xlabel('$T$', 'interpreter', 'latex')
ylabel('$L(q)$', 'interpreter', 'latex');
print('varL_m1', '-depsc2', '-r300');
% eps2xxx('varL_1.eps', {'jpeg'}, 'C:\Program Files\gs\gs9.10\bin\gswin64c.exe')