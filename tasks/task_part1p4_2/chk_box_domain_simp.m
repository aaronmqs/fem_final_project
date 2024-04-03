function chk_box_domain_simp(ndim)
% This function checks if the volume, the centroid and the surface area of
% a ndim-dimensional domain which mesh is formed by simplicial elements is
% evaluated correctly.

% Create mesh
etype = 'simp';
porder = 1;
lims = [zeros(ndim, 1) ones(ndim, 1)];
nel = 2^ndim * ones(ndim, 1);
msh_quad = create_mesh_hcube('hcube', lims, nel, porder);
[e2vcg_tri, e2bnd_tri] = split_quad_mesh_into_tri_mesh(msh_quad.e2vcg, msh_quad.e2bnd);
msh = create_mesh_strct(etype, msh_quad.xcg, e2vcg_tri, e2bnd_tri);

% Visualize mesh
if ndim == 1 || ndim == 2
    visualize_fem([], msh);
elseif ndim == 3
    visualize_fem3d([], msh);
else
    warning("Visualization not implemented.")
end

% Compute metrics (volume, centroid, surface area)
nquad_per_dim = ceil((porder + 1) / 2) * 3; % TODO: Why so many qpoints necessary?
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
fprintf("Test for simplex passed.\n\n")

end