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

% Preallocate element residual and Jacobian
Re = zeros(nvar_per_elem, 1);
dRe = zeros(nvar_per_elem, nvar_per_elem);

% Code me!
Tv_ref = elem.Tv_ref;
for q = 1:nq
    Phie = Tv_ref(:, :, 1, q);
    U(:, q) = Phie' * Ue;
    for j = 1:ndim
        Q(:, j, q) = Tv_ref(:, :, j, q)' * Ue;
    end
end

for l = 1:nvar_per_elem
    for q = 1:nq
        Psi = Tv_ref(l, :, 1, q);
        dPsi = Tv_ref(l, :, 2:end, q);
%         Re(l, 1) = Re(l, 1) +  () * detG(nq) * wq(nq);
    end
end























end