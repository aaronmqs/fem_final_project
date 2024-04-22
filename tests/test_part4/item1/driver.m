clc; clear; close all;

% Create equation
ndim = 2;
eqn = IncompressibleNavierStokes(ndim);
vol_pars_fcn = @(x) [1; 0.1];
bnd_pars_fcn = @(x, bnd) [-1 -2 0]';

% Extract information from input
nvar = eqn.nvar;

% Create a scalar local function space
etype = 'simp'; nquad_per_dim = 10; porder = 2;
qrule = create_qrule_gaussleg(etype, ndim, nquad_per_dim);
lfcnsp = create_polysp_nodal(etype, ndim, porder, qrule.zq, qrule.rq);

% Create element basis functions
Qv = lfcnsp.Qv; Qvf = lfcnsp.Qvf;
[Tv_ref, Tvf_ref] = create_elem_basis(nvar, Qv, Qvf);

% Create ELEM structure
elem = create_elem_strct(eqn, qrule, lfcnsp);

% Create mesh
c = zeros(ndim, 1); r = 1; nel = 5 * ones(ndim, 1);
msh = create_mesh_hsphere(etype, c, r, nel, porder);
% visualize_fem([], msh);

% Compute transformation quantities
transf_data = create_transf_data_ndim(lfcnsp, msh.xcg, msh.e2vcg, msh.e2bnd);

% Create elem_data structure
elem_data = create_elem_data(Tv_ref, Tvf_ref, ...
                                        transf_data, vol_pars_fcn, ...
                                        bnd_pars_fcn);

% Create random state
ndof_per_elem = size(Tv_ref, 1);
Ue = rand(ndof_per_elem, 1);

% Integrate element Galerkin form (volume term) and Jacobian (element 1)
[Re_vol, dRe] = intg_elem_claw_vol(Ue, transf_data(1), elem(1), elem_data(1));

Re_bnd = intg_elem_claw_extface(transf_data(1), elem(1), elem_data(1));

Re = Re_vol + Re_bnd;


