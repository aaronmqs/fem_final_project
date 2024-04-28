function [U, info] = solve_linelptc_sclr_nd0(nref, porder, pltit)
%SOLVE_LINELPTC_SCLR_ND0 Solve Poisson equation on ND mesh using
%FEM with simplex elements of polynomial completeness PORDER, see project
%handout for complete problem description.
%
% Input arguments
% ---------------
%   NREF : number : level of refinement (only 0 supported for now)
%
%   PORDER : See notation.m
%
%   PLTIT : bool : Whether to plot solution
%
% Output arguments
% ----------------
%   U : Array (NDOF,) : Global (assembled) finite element solution
%
%   INFO : See NEWTRAPH

ndim = 2;
nvar = 1;

% Load finite element mesh
msh = load_mesh('nd0', 'simp', nref, porder);
xcg = msh.xcg; e2vcg = msh.e2vcg; e2bnd = msh.e2bnd;

% Setup equation parameters and natural boundary conditions
prob.eqn = LinearEllipticScalar(2);
% prob.vol_pars_fcn = % TODO
prob.vol_pars_fcn = @(x) [1; 0; 0; 1; 0];
% prob.bnd_pars_fcn = % TODO
prob.bnd_pars_fcn = @(x, bnd) 0;

% Extract indices and set values of dirichlet boundary conditions
[~, f2v, ~] = create_nodes_bndy_refdom_simp(ndim, porder);

% dbc_idx = % TODO
ndof_per_node = nvar;
ldof2gdof = create_ldof2gdof_cg(ndof_per_node, e2vcg);
ldof = 1;
dbc_idx1 = get_gdof_from_bndtag(ldof, 1, ndof_per_node, ldof2gdof, e2bnd, f2v);
dbc_idx2 = get_gdof_from_bndtag(ldof, 2, ndof_per_node, ldof2gdof, e2bnd, f2v);
dbc_idx = get_gdof_from_bndtag(ldof, [1 2], ndof_per_node, ldof2gdof, e2bnd, f2v);
% dbc_val = % TODO
dbc_val = [zeros(size(dbc_idx1)); 10 * ones(size(dbc_idx2))];

ndof = size(xcg, 2);
dbc = create_dbc_strct(ndof, dbc_idx, dbc_val);

% Create finite element space
femsp = create_femsp_cg(prob, msh, porder, e2vcg, dbc);

% Solve finite element equations
tol = 1.0e-8; maxit = 10;
[U, info] = solve_fem(femsp, [], tol, maxit);
               
% Visualize FEM solution
if pltit
    % Evaluate FEM solution throughout domain
%     xeval = % TODO
    xeval1 = [linspace(-1/6, 7/6); -0.25 * ones(1, 100)];
    Ux1 = eval_fem_soln(U(femsp.ldof2gdof), xeval1, msh, femsp.elem);

    xeval2 = [linspace(-1/6, 7/6); -0.64 * ones(1, 100)];
    Ux2 = eval_fem_soln(U(femsp.ldof2gdof), xeval2, msh, femsp.elem);

    xeval3 = [0.33 * ones(1, 100); linspace(-9/8, 1/8)];
    Ux3 = eval_fem_soln(U(femsp.ldof2gdof), xeval3, msh, femsp.elem);

    xeval4 = [0.73 * ones(1, 100); linspace(-9/8, 1/8)];
    Ux4 = eval_fem_soln(U(femsp.ldof2gdof), xeval4, msh, femsp.elem);

    visualize_fem([], msh, U(e2vcg), struct('plot_elem', true, 'nref', 2));
    colorbar;

    figure;
    plot(xeval1(1, :), squeeze(Ux1(1, 1, :)), 'k-', 'linewidth', 2);

    figure
    plot(xeval2(1, :), squeeze(Ux2(1, 1, :)), 'k-', 'linewidth', 2);

    figure
    plot(xeval3(1, :), squeeze(Ux3(1, 1, :)), 'k-', 'linewidth', 2);

    figure
    plot(xeval4(1, :), squeeze(Ux4(1, 1, :)), 'k-', 'linewidth', 2);
end

end