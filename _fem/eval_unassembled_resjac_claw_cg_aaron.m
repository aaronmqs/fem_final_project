function  [Re, dRe] = eval_unassembled_resjac_claw_cg(U, transf_data, elem, elem_data, ldof2gdof)
%EVAL_UNASSEMBLED_RESJAC_CLAW_CG Evaluate/store element residual vector and
%Jacobian matrices for each element. 
%
% Input arguments
% ---------------
%   U : Array (NDOF,) : Global (assembled) solution vector
%
%   TRANSF_DATA, ELEM, ELEM_DATA, LDOF2GDOF : See notation.m 
%
% Output arguments
% ----------------
%   RE : Array (NDOF_PER_ELEM, NELEM): Element residual vector for
%     all elements in mesh
%
%   DRE : Array (NDOF_PER_ELEM, NDOF_PER_ELEM, NELEM): Element Jacobian
%     matrix for all elements in mesh

% Extract relevant variables from element
nelem = numel(elem_data);
nvar_per_elem = size(elem.Tv_ref, 1);
if isempty(U), ndof = max(ldof2gdof(:)); U = zeros(ndof, 1); end
    
% Preallocate and assemble
Re = zeros(nvar_per_elem, nelem);
dRe = zeros(nvar_per_elem, nvar_per_elem, nelem);

% Code me!
for e = 1:nelem
    idxe = ldof2gdof(:, e);
    Ue = U(idxe);
    [Re_vol, dRe(:, :, e)] = intg_elem_claw_vol(Ue, transf_data(e), elem, elem_data(e));
    Re_bnd = intg_elem_claw_extface(transf_data(e), elem, elem_data(e));
    Re(:, e) = Re_vol + Re_bnd;
end

end