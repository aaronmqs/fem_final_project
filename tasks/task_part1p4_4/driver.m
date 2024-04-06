clear; clc; close all;

porder = 1; etype = 'simp';

% 2d meshes
ndim2 = 2; 

% Batman domain
msh_batman = load_mesh('batman0', etype, 0, porder);
[vb, cb, sab] = measure_domain(msh_batman, porder, ndim2, etype);
fname = sprintf('C:/Users/aaron/PhD/2024_1_AME60541/project/fem_final_project/tasks/task_part1p4_4/batman');
saveas(gcf, [fname, '.png']); close all;

% ND domain
msh_nd = load_mesh('nd0', etype, 0, porder);
[vnd, cnd, sand] = measure_domain(msh_nd, porder, ndim2, etype);
fname = sprintf('C:/Users/aaron/PhD/2024_1_AME60541/project/fem_final_project/tasks/task_part1p4_4/nd');
saveas(gcf, [fname, '.png']); close all;

% 3d meshes
ndim3 = 3;

% Cow domain
msh_cow = load_mesh('cow', etype, 0, porder);
[vc, cc, sac] = measure_domain(msh_cow, porder, ndim3, etype);
fname = sprintf('C:/Users/aaron/PhD/2024_1_AME60541/project/fem_final_project/tasks/task_part1p4_4/cow');
saveas(gcf, [fname, '.png']); close all;

% Dragon domain
msh_dragon = load_mesh('dragon', etype, 0, porder);
[vd, cd, sad] = measure_domain(msh_dragon, porder, ndim3, etype);
fname = sprintf('C:/Users/aaron/PhD/2024_1_AME60541/project/fem_final_project/tasks/task_part1p4_4/dragon');
saveas(gcf, [fname, '.png']); close all;

% Skultp10kv domain
msh_skultp10kv = load_mesh('sculpt10kv', etype, 0, porder);
[vs, cs, sas] = measure_domain(msh_skultp10kv, porder, ndim3, etype);
fname = sprintf('C:/Users/aaron/PhD/2024_1_AME60541/project/fem_final_project/tasks/task_part1p4_4/sculpt10kv');
saveas(gcf, [fname, '.png']); close all;