function [Ru, dRu] = create_fem_resjac(Uu, femsp)
%CREATE_FEM_RESJAC Create the finite element residual and Jacobian,
%restricted to the free degrees of freedom. When combined with a nonlinear
%solver, this will approximate the solution of a PDE (described in FEMSP).
%
% Input arguments
% ---------------
%   UU : Array (NDOF-NDBC,) : Global (assembled) solution vector,
%     restricted to the free degrees of freedom (via static condensation).
%
%   FEMSP : See notation.m
%
% Output arguments
% -----------------
%   RU : Array (NDOF-NDBC,) : Finite element residual, restricted to free
%     degrees of freedom
%
%   DRU : Sparse matrix (NDOF-NDBC, NDOF-NDBC) : Finite element Jacobian
%     restricted to free degrees of freedom

% Extract information from input
ndof = max(femsp.ldof2gdof(:));

% Code me!
dbc = femsp.dbc;
dbc_idx = dbc.dbc_idx;
free_idx = dbc.free_idx;
dbc_val = dbc.dbc_val;
U(dbc_idx, 1) = dbc_val;
U(free_idx, 1) = Uu;
[R, dR] = eval_assembled_resjac_claw_cg(U, femsp.transf_data, femsp.elem, femsp.elem_data, ...
                                                 femsp.ldof2gdof, femsp.spmat);

Ru = R(free_idx, 1);
dRu = dR(free_idx, free_idx);
end