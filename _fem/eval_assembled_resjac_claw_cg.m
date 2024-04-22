function [R, dR] = eval_assembled_resjac_claw_cg(U, transf_data, elem, elem_data, ...
                                                 ldof2gdof, spmat)
%EVAL_ASSEMBLED_RESJAC_CLAW_CG Evaluate assembled residual vector and
%Jacobian matrix
%
%Input arguments
%---------------
%   U : Array (NDOF,) : Global (assembled) solution vector
%
%   TRANSF_DATA, ELEM, ELEM_DATA, SPMAT : See notation.m
%
%Output arguments
%----------------
%   R : Array (NDOF,) : Assembled residual vector PRIOR to static condensation
%
%   dR : Array (NDOF, NDOF) : Assembled Jacobian matrix PRIOR to static condensation

% Code me!
cooidx = spmat.cooidx;
lmat2gmat = spmat.lmat2gmat;

[Ru, dRu] = eval_unassembled_resjac_claw_cg(U, transf_data, elem, elem_data, ldof2gdof);
dR = assemble_nobc_mat(dRu, cooidx, lmat2gmat);
R = assemble_nobc_vec(Ru, ldof2gdof);

end