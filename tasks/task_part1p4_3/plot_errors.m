clear; clc; close all;

% Load errors
load("ec.mat"); load("esa.mat"); load("ev.mat");

etype = {'hcube', 'simp'};
for ndim = 2:3
    ndimm1 = ndim - 1;
    nelem_per_dim = 1:2:2^(3 + 1);
    for i = 1:2 % 1=hcube, 2=simp
        if ~(i == 2 && ndim == 3)
            figure
            for porder = 4:4
                idx = 1:length(nelem_per_dim);
                subplot(3, 1, 1)
                hold on;
                loglog(nelem_per_dim(idx), squeeze(ec(i, ndimm1, porder, idx)), 'Marker','o'); grid on; 
                ylabel("Centroid error")
                xticks(nelem_per_dim)
                subplot(3, 1, 2)
                hold on;
                loglog(nelem_per_dim(idx), squeeze(esa(i, ndimm1, porder, idx)), 'Marker', 'diamond'); grid on;
                ylabel("Surf. area error")
                xticks(nelem_per_dim)
                subplot(3, 1, 3)
                hold on;
                semilogy(nelem_per_dim(idx), squeeze(ev(i, ndimm1, porder, idx)), 'Marker', 'square'); grid on;
                ylabel("Volume error")
                xlabel("Number of elements");
                xticks(nelem_per_dim)
            end
            str = sprintf("%d-dimensional %s element (p = 4).", ndim, etype{i});
            sgtitle(str);
            fname = sprintf('C:/Users/aaron/PhD/2024_1_AME60541/project/fem_final_project/tasks/task_part1p4_3/nel_ndim%g_%s', ndim, etype{i});
            saveas(gcf, [fname, '.png']); close all;
        end
    end
end

for ndim = 2:3
    ndimm1 = ndim - 1;
    nelem_per_dim = 1:2:2^(3 + 1);
    for i = 1:2 % 1=hcube, 2=simp
        if ~(i == 2 && ndim == 3)
            figure
            for idx = length(nelem_per_dim):length(nelem_per_dim)
                subplot(3, 1, 1)
                hold on;
                loglog(1:4, squeeze(ec(i, ndimm1, :, idx)), 'Marker','o'); grid on; 
                ylabel("Centroid error");
                xticks([1 2 3 4])
                subplot(3, 1, 2)
                hold on;
                loglog(1:4, squeeze(esa(i, ndimm1, :, idx)), 'Marker', 'diamond'); grid on;
                ylabel("Surf. area error");
                xticks([1 2 3 4])
                subplot(3, 1, 3)
                hold on;
                semilogy(1:4, squeeze(ev(i, ndimm1, :, idx)), 'Marker', 'square'); grid on;
                ylabel("Volume error")
                xlabel("Polynomial degree");
                xticks([1 2 3 4])
            end
            str = sprintf("%d-dimensional %s element (nelem = %d).", ndim, etype{i}, nelem_per_dim(end));
            sgtitle(str);
            fname = sprintf('C:/Users/aaron/PhD/2024_1_AME60541/project/fem_final_project/tasks/task_part1p4_3/p_ndim%g_%s', ndim, etype{i});
            saveas(gcf, [fname, '.png']); close all;
        end
    end
end

