clear; clc; close all;

% diary log.txt

bffr = NaN(2, 2, 4, 8);
[ev, esa, ec] = deal(bffr, bffr, bffr);
i = 0; tol = 1e-6; pltit = 0;
for etype = {'hcube', 'simp'}
    i = i + 1;
    for ndim = 2:3
        ndimm1 = ndim - 1;
        for porder = 1:4            
            for j = 1:2^3 % more points because of the mappings (same number for any dimensions)
                nel = 2 * j * ones(ndim, 1);
                if ~(strcmp(etype{:}, 'simp') && ndim == 3)
                    fprintf("\nTest for %d-dimensional mesh of %d %s elements (per dimension) of order %d: \n\n", ndim, i * nel(1), etype{:}, porder);
                    [ev(i, ndimm1, porder, j), esa(i, ndimm1, porder, j), ec(i, ndimm1, porder, j)] = chk_sphere_domain(etype{:}, ndim, nel, porder, pltit);
                    fprintf("Volume error: %d. \n", ev(i, ndimm1, porder, j))
                    fprintf("Surface area error: %d. \n", esa(i, ndimm1, porder, j))
                    fprintf("Centroid error: %d. \n", ec(i, ndimm1, porder, j))
                    if norm([ev(i, ndimm1, porder, j), esa(i, ndimm1, porder, j), ec(i, ndimm1, porder, j)], "inf") < tol
                        fprintf("\nThis number of elements/polynomial order are enough to evaluate the integrals for this mesh. \n")
                    end
                end
            end
        end
    end
end

% save("ev", "ev")
% save("esa", "esa")
% save("ec", "ec")
% 
% diary off