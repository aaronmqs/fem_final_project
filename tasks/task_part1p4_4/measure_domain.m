function [v, c, sa] = measure_domain(msh, porder, ndim, etype)
% This function computes the volume, the centroid and the surface area of
% a mesh defined domain.

% Compute metrics (volume, centroid, surface area)
nquad_per_dim = ceil((porder + 1) / 2) * 3; % More quadrature points than "usual" will be necessary, since the elements were transformed.
qrule = create_qrule_gaussleg(etype, ndim, nquad_per_dim);
lfcnsp = create_polysp_nodal(etype, ndim, porder, qrule.zq, qrule.rq);
transf_data = create_transf_data_ndim(lfcnsp, msh.xcg, msh.e2vcg, msh.e2bnd);
[v, c, sa] = compute_domain_metrics(transf_data, qrule);

% % Visualize mesh
% visualize_fem([], msh);
% if ndim == 2
%     scatter(c(1), c(2), 'filled', 'MarkerFaceColor','k')
% elseif ndim == 3
%     scatter3(c(1), c(2), c(3), 'filled', 'MarkerFaceColor','k')
% end

end