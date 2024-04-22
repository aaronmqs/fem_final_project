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

% Preallocate
S = zeros(neqn, 1, nq);
dSdU = zeros(neqn, nvar, nq);
dSdQ = zeros(neqn, nvar, ndim, nq);
F = zeros(neqn, ndim, nq);
dFdU = zeros(neqn, ndim, nvar, nq);
dFdQ = zeros(neqn, ndim, nvar, ndim, nq);
U = zeros(nvar, nq);
Q = zeros(nvar, ndim, nq);

% Solution basis evaluated at the quadrature nodes
Tv_ref = elem.Tv_ref;
Tv_phys = elem_data.Tv_phys;
vol_pars = elem_data.vol_pars;

for q = 1:nq
    % State gradient
    Phie = Tv_ref(:, :, 1, q);
    for j = 1:ndim
        Q(:, j, q) = Tv_phys(:, :, j, q)' * Ue;
    end

    % Source and flux
    [S(:, q), dSdU(:, :, q), dSdQ(:, :, :, q), F(:, :, q), dFdU(:, :, :, q), dFdQ(:, :, :, :, q)] = elem.eqn.srcflux(Phie' * Ue, Q(:, :, q), vol_pars(:, q));

    Psi = Tv_ref(:, :, 1, q);
    dPsi = squeeze(Tv_phys(:, :, 2:end, q));

    part1 = dSdU(:, :, q) * Psi';
    part2 = zeros(size(part1));
    part3 = zeros(neqn, ndim, nvar_per_elem);
    part4 = zeros(size(part3));
    part6 = zeros(nvar_per_elem, nvar_per_elem);
    for r = 1:nvar_per_elem
        for i = 1:neqn
            part2(i, r) = sum( squeeze(dSdQ(i, :, :, q)) .* squeeze(dPsi(r, :, :) ) , 'all');
            part3(i, :, r) = squeeze(dFdU(i, :, :, q)) * Psi(r, :)';
            for j = 1:ndim
                part4(i, j, r) = sum(squeeze(dFdQ(i, j, :, :, q)) .* squeeze(dPsi(r, :, :)), 'all');
            end
        end
        part5 = part3 + part4;
        for l = 1:nvar_per_elem
            part6(l, r) = - sum(squeeze(dPsi(l, :, :)) .* part5, 'all');
        end
        Re(r, 1) = Re(r, 1) + (- Psi(r, :) * S(:, q) - sum(squeeze(dPsi(r, :, :)) .* F(:, :, q), 'all')) * detG(q) * wq(q);
    end
    part7 = part1 + part2;
    part8 = - Psi * part7;
    integrand = (part8 + part6) * detG(q);
    dRe = dRe + integrand * wq(q);
end

end