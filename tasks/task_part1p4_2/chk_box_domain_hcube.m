function chk_box_domain_hcube(ndim)
% This function checks if the volume, the centroid and the surface area of
% a ndim-dimensional domain which mesh is formed by hypercube elements is
% evaluated correctly.

% Create mesh
etype = 'hcube';
porder = 1;
lims = [zeros(ndim, 1) ones(ndim, 1)];
nel = 2^ndim * ones(ndim, 1);
msh = create_mesh_hcube(etype, lims, nel, porder);

% Visualize mesh
if ndim == 1 || ndim == 2
    visualize_fem([], msh);
elseif ndim == 3
    visualize_fem3d([], msh);
else
    warning("Visualization not implemented.")
end

% Compute metrics (volume, centroid, surface area)
nquad_per_dim = ceil((porder + 1) / 2); % Only this number of qpoints is needed, since the elements are all perfect hypercubes.
qrule = create_qrule_gaussleg(etype, ndim, nquad_per_dim);
lfcnsp = create_polysp_nodal(etype, ndim, porder, qrule.zq, qrule.rq);
transf_data = create_transf_data_ndim(lfcnsp, msh.xcg, msh.e2vcg, msh.e2bnd);
[v, c, sa] = compute_domain_metrics(transf_data, qrule);

% Checks
tol = 1e-8;
if abs(v - 1) > tol; error("Incorrect volume."); end
if ~isempty(find(abs(c - 0.5) > tol, 1)); error("Incorrect centroid."); end
if abs(sa - 2 * ndim) > tol; error("Incorrect surface area."); end

% If this point is reached, the test was successful.
fprintf("Test for hypercube passed.\n\n")

end