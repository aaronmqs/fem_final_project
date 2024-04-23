clear; clc; close all;

% Create problem
prob.eqn = Pde0();
prob.vol_pars_fcn = @(x) x^2;
prob.bnd_pars_fcn = @(x, bnd) -1 * (x >= 1 - 1e-12);

% Create mesh
etype = 'hcube'; lims = [0 5]; nel = 50; porder = 2;
msh = create_mesh_hcube(etype, lims, nel, porder);
% visualize_fem([], msh);

% Create Dirichlet boundary conditions
ndof = numel(msh.xcg); dbc_idx = 1; dbc_val = 0.0;
dbc = create_dbc_strct(ndof, dbc_idx, dbc_val);

% Create finite element space
e2vcg_porder = msh.e2vcg;
femsp = create_femsp_cg(prob, msh, porder, e2vcg_porder, dbc);

% Test function that creates residual and its jacobian retricted to the
% unconstrained DoFs.
Uu = rand();
[Ru, dRu] = create_fem_resjac(Uu, femsp);