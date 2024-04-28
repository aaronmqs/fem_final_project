function [Re, dRe] = intg_elem_claw_vol(Ue, transf_data, elem, elem_data)
%INTG_ELEM_CLAW_VOL Integrate element Galerkin form (volume term) to
%form the volume contribution to the element residual and Jacobian.
%
% Input arguments
% ---------------
%   UE : Array (NDOF_PER_ELEM,) : Element solution (primary variables)
%
%   TRANSF_DATA, ELEM, ELEM_DATA : See notation.m
%
% Output arguments
% ----------------
%   RE : Array (NDOF_PER_ELEM,) : Element residual (volume contribution)
%
%   DRE : Array (NDOF_PER_ELEM, NDOF_PER_ELEM) : Element Jacobian (volume contribution)

% Extract information from input : sizes
[nvar_per_elem, nvar, ndimP1, nq] = size(elem.Tv_ref);
neqn = nvar; ndim = ndimP1-1;

% Extract information from input : quadrature
wq = elem.qrule.wq;

% Extract information from input : isoparametric
detG = transf_data.detG;

% Extract information from input : basis
Tvar = reshape(elem_data.Tv_phys, [nvar_per_elem, nvar*(ndim+1)*nq]);

% Preallocate element residual and Jacobian
Re = zeros(nvar_per_elem, 1);
dRe = zeros(nvar_per_elem, nvar_per_elem);

% Move primary variables from nodes to quadrature nodes
UQq = reshape(Tvar'*Ue, [nvar, ndim+1, nq]);

% Integrate volume term
w = wq.*detG;
for k = 1:nq
    % Solution basis functions and a derivative at quadrature point
    Tvar = reshape(elem_data.Tv_phys(:, :, :, k), [nvar_per_elem, nvar*(ndim+1)]);
    
    % Evaluate pointwise quantities (flux, source) and reshape
    pars = elem_data.vol_pars(:, k);
    [S, dSdU, dSdQ, F, dFdU, dFdQ] = elem.eqn.srcflux(UQq(:, 1, k), UQq(:, 2:end, k), pars);
    SF = [S, F]; SF = SF(:);
    dSFdU = cat(2, reshape(dSdU, [neqn, 1, nvar]), dFdU);
    dSFdQ = cat(2, reshape(dSdQ, [neqn, 1, nvar, ndim]), dFdQ);
    dSFdUQ = cat(4, reshape(dSFdU, [neqn, ndim+1, nvar, 1]), dSFdQ);
    dSFdUQ = reshape(dSFdUQ, [neqn*(ndim+1), nvar*(ndim+1)]);

    % Add contribution to element residual and Jacobian
    Re = Re - w(k)*(Tvar*SF);
    dRe = dRe - w(k)*(Tvar*(dSFdUQ*Tvar'));
end

end