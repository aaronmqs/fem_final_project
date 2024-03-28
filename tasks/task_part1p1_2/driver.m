clear; clc; close all;

d = 2; p = 2;
[zk, f2v, N] = create_nodes_bndy_refdom_simp(d, p);

nf = size(f2v, 2);
mrkrs = {'x'; 'square'; '+'};
clrs = ['r', 'g', 'k'];
for i = 1:nf
    figure
    plot(zeros(100, 1), linspace(0, 1, 100), 'k')
    hold on
    plot(linspace(0, 1, 100), zeros(100, 1), 'k')
    hold on
    plot(linspace(0, 1, 1000), 1 - linspace(0, 1, 1000), 'k')
    hold on
    scatter(zk(1,:), zk(2,:), 'filled', 'MarkerFaceColor', 'b')
    hold on
    nds_fc = f2v(:, i);
    scatter(zk(1, nds_fc), zk(2, nds_fc), clrs(i), 'Marker', mrkrs{i}, 'LineWidth', 2, 'SizeData', 300)
    title('Nodes on face', i)
    grid on
    axis([-0.1 1.1 -0.1 1.1])

%     fname = sprintf('C:/Users/aaron/PhD/2024_1_AME60541/project/fem_final_project/tasks/Part1p1/face%d', i); saveas(gcf, [fname, '.png']); close all;
end