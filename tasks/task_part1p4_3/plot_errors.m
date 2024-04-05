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
                idx = 4:length(nelem_per_dim);
                subplot(3, 1, 1)
                hold on;
                semilogy(nelem_per_dim(idx), squeeze(ec(i, ndimm1, porder, idx))); grid on;
                ylabel("Centroid error")
                legend("p1", "p2", "p3", "p4")
                subplot(3, 1, 2)
                hold on;
                semilogy(nelem_per_dim(idx), squeeze(esa(i, ndimm1, porder, idx))); grid on;
                ylabel("Surf. area error")
                legend("p1", "p2", "p3", "p4")
                subplot(3, 1, 3)
                hold on;
                semilogy(nelem_per_dim(idx), squeeze(ev(i, ndimm1, porder, idx))); grid on;
                ylabel("Volume error")
                legend("p1", "p2", "p3", "p4")
                xlabel("Number of elements");
            end
            str = sprintf("%d-dimensional %s element.", ndim, etype{i});
            sgtitle(str);
        end
    end
end