vS = 20;
hS = 20;
T = 0.5 : 0.1 : 10;
betaAll = 1 ./ T;
H = zeros(vS, hS);
J = -1;
num_iter = 10000;
opt_params.max_iter = 100;
opt_params.num_start = 20;

tic
[E4, D4, M4, S4] = gibbsIsing(H, J, betaAll, num_iter, 4);
[E6, D6, M6, S6] = gibbsIsing(H, J, betaAll, num_iter, 6);
toc

figure()
hold on
plot(T, D4, 'b')
plot(T, D6, 'r')
legend('rectangular', 'triangular')
xlabel('$T$', 'interpreter', 'latex')
ylabel('$\frac{1}{N}\sqrt{\mathbf{D}[X]}$', 'interpreter', 'latex');%, 'FontSize', 13
print('D_m1', '-depsc2', '-r300');
% eps2xxx('D_m1.eps', {'jpeg'}, 'C:\Program Files\gs\gs9.10\bin\gswin64c.exe')

figure()
hold on
plot(T, E4, 'b')
plot(T, E6, 'r')
legend('rectangular', 'triangular', 'Location', 'SouthEast')
xlabel('$T$', 'interpreter', 'latex')
ylabel('$\frac{1}{N}\mathbf{E}[X]$', 'interpreter', 'latex');
print('E_m1', '-depsc2', '-r300');
% eps2xxx('E_m1.eps', {'jpeg'}, 'C:\Program Files\gs\gs9.10\bin\gswin64c.exe')

figure()
hold on
plot(T, M4, 'b')
plot(T, M6, 'r')
legend('rectangular', 'triangular', 'Location', 'SouthEast')
xlabel('$T$', 'interpreter', 'latex')
ylabel('$\sqrt{\mathbf{E}[\mu^2(X)]}$', 'interpreter', 'latex');
print('M_m1', '-depsc2', '-r300');
% eps2xxx('M_m1.eps', {'jpeg'}, 'C:\Program Files\gs\gs9.10\bin\gswin64c.exe')