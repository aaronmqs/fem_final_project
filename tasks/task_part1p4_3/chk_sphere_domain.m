function [varargout] = chk_sphere_domain(etype, ndim, nel, porder)
% This function checks if the volume, the centroid and the surface area of
% a ndim-dimensional sphere domain which mesh is formed by hypercube elements is
% evaluated correctly.

if ndim ~= 2 && ndim ~= 3; error("Not implmented."); end

% Create mesh
c = zeros(ndim, 1); r = 1;
msh = create_mesh_hsphere(etype, c, r, nel, porder);

% Visualize mesh
if ndim == 1 || ndim == 2
    visualize_fem([], msh);
elseif ndim == 3
    visualize_fem3d([], msh);
else
    warning("Visualization not implemented.")
end

% Compute metrics (volume, centroid, surface area)
nquad_per_dim = ceil((porder + 1) / 2) * 3; % More quadrature points than "usual" will be necessary, since the elements were transformed.
qrule = create_qrule_gaussleg(etype, ndim, nquad_per_dim);
lfcnsp = create_polysp_nodal(etype, ndim, porder, qrule.zq, qrule.rq);
transf_data = create_transf_data_ndim(lfcnsp, msh.xcg, msh.e2vcg, msh.e2bnd);
[v, c, sa] = compute_domain_metrics(transf_data, qrule);

% Checks
if nargout > 0
    if ndim == 2
        ev = abs(v - pi * r^2);
        esa = abs(sa - 4 * pi * r^2);
    else
        ev = abs(v - 4 * pi * r^3 / 3);
        esa = abs(sa - 4 * pi * r^2);
    end
    ec = norm(c);
    varargout = [ev, esa, ec];
else
    tol = 1e-8;
    if ndim == 2
        if abs(v - pi * r^2) > tol; error("Incorrect volume."); end
        if abs(sa - 2 * pi * r) > tol; error("Incorrect surface area."); end
    else
        if abs(v - 4 * pi * r^3 / 3) > tol; error("Incorrect volume."); end
        if abs(sa - 4 * pi * r^2) > tol; error("Incorrect surface area."); end
    end    

    if ~isempty(find(abs(c) > tol, 1)); error("Incorrect centroid."); end
end


% If this point is reached, the test was successful.
fprintf("Test for hypercube passed.\n\n")

end