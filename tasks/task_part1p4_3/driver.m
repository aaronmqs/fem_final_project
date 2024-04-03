clear; clc; close all;

% hcube
for ndim = 2:3
    for porder = 1:4
        chk_sphere_domain_hcube(ndim, nel, porder);
    end
end
