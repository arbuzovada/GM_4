% S4(S4 == -1) = 0;
% for i = 1 : 96
%     imwrite(repmat(S4(:,:,i), [1,1,3]), ['.\report\figures\S_', num2str(i), '.png'])
% end

T = 0.5 : 0.1 : 10;
betaAll = 1 ./ T;

f = figure('Position', [100, 100, 500, 500]);
load q_final.mat
q = best_q;
q = reshape(q, [vS, hS, beta0]);
for i = 1 : beta0
    image(repmat(S4(:, :, i), [1, 1, 3]))
    axis equal
    box on
    set(gca, 'XLim', [0.5, hS + 0.5]);
    set(gca, 'YLim', [0.5, vS + 0.5]);
    set(gca, 'XTick', [])
    set(gca, 'YTick', [])
    title(['$T$ = ', num2str(T(i))], 'interpreter', 'latex')
    saveas(f, ['.\report\figures\S_', num2str(i), '.png'])
end
% print('q', '-depsc2', '-r300');
% eps2xxx('q.eps', {'png'}, 'C:\Program Files\gs\gs9.10\bin\gswin64c.exe')