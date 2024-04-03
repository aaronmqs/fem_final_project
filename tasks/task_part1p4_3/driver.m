clear; clc; close all;

bffr = zeros(2, 2, 4);
[ev, esa, ec] = deal(bffr, bffr, bffr);
i = 0;
tol = 1e-6; pltit = 0;
for etype = {'hcube', 'simp'}
    i = i + 1;
    for ndim = 2:3
        for porder = 1:4
            nel_base = 2^ndim * ones(ndim, 1);
            for j = 1:2
                nel = j * nel_base;
                if ~(strcmp(etype{:}, 'simp') && ndim == 3)
                    fprintf("\nTest for %d-dimensional mesh of %d %s elements (per dimension) of order %d: \n\n", ndim, i * nel(1), etype{:}, porder);
                    [ev(i, ndim, porder, j), esa(i, ndim, porder, j), ec(i, ndim, porder, j)] = chk_sphere_domain(etype{:}, ndim, nel, porder, pltit);
                    fprintf("Volume error: %d. \n", ev(i, ndim, porder, j))
                    fprintf("Surface area error: %d. \n", esa(i, ndim, porder, j))
                    fprintf("Centroid error: %d. \n", ec(i, ndim, porder, j))
                    if norm([ev(i, ndim, porder, j), esa(i, ndim, porder, j), ec(i, ndim, porder, j)], "inf") < tol
                        fprintf("\nThis number of elements/polynomial order are enough to evaluate the integrals for this mesh. \n")
                    end
                end
            end
        end
    end
end

save("ev", "ev")
save("esa", "esa")
save("ec", "ec")