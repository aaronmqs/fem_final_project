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
Tv_phys = elem_data.Tv_phys;

for q = 1:nq
    U = Tv_phys(:, :, 1, q)' * Ue;
    dUdX = zeros(neqn, ndim);
    for j = 1:ndim
        dUdX(:, j) = Tv_phys(:, :, 1 + j, q)' * Ue;
    end
    pars = elem_data.vol_pars(:, q);
    [S, dSdU, dSdQ, F, dFdU, dFdQ] = elem.eqn.srcflux(U, dUdX, pars);
    
    % Residual
    integrand_Re = zeros(nvar_per_elem, 1);
    for l = 1:nvar_per_elem
        integrand_Re(l, 1) = - Tv_phys(l, :, 1, q) * S;
        for j = 1:ndim
            integrand_Re(l, 1) = integrand_Re(l, 1) - Tv_phys(l, :, 1 + j, q) * F(:, j);
        end
    end
    integrand_Re = integrand_Re * detG(q);
    Re = Re + integrand_Re * wq(q);
    
    % Jacobian
    integrand_dRe = zeros(nvar_per_elem, nvar_per_elem);
    for l = 1:nvar_per_elem
        for r = 1:nvar_per_elem
            integrand_dRe(l, r) = - Tv_phys(l, :, 1, q) * dSdU * Tv_phys(r, :, 1, q)';
            for s = 1:ndim
                integrand_dRe(l, r) = integrand_dRe(l, r) - Tv_phys(l, :, 1, q) * dSdQ(:, :, s) * Tv_phys(r, :, 1 + s, q)';
            end
            for j = 1:ndim
                integrand_dRe(l, r) = integrand_dRe(l, r) - Tv_phys(l, :, 1 + j, q) * squeeze(dFdU(:, j, :)) * Tv_phys(r, :, 1, q)';
                for s = 1:ndim
                    integrand_dRe(l, r) = integrand_dRe(l, r) - Tv_phys(l, :, 1 + j, q) * squeeze(dFdQ(:, j, :, s)) * Tv_phys(r, :, 1 + s, q)';
                end
            end
        end
    end
    integrand_dRe = integrand_dRe * detG(q);
    dRe = dRe + integrand_dRe * wq(q);
end

end