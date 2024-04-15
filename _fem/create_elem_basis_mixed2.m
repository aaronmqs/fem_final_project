function [Tv, Tvf] = create_elem_basis_mixed2(nvar1, Qv1, Qvf1, ...
                                              nvar2, Qv2, Qvf2)
%CREATE_ELEM_BASIS_MIXED2 Create element basis for mixed element (two
%different local function spaces) evaluated at point throughout volume (TV)
%and on each face (TVF) from basis of local function spaces evaluated at
%corresponding points (QV1, QVF1, QV2, QVF2). 
%
% Input arguments
% ---------------
%   NVAR1, QV1, QVF1 : NVAR, QV, QVF (see notation.m) for the first local
%     function space, e.g., velocity in Navier-Stokes.
%
%   NVAR2, QV2, QVF2 : NVAR, QV, QVF (see notation.m) for the second local
%     function space, e.g., pressure in Navier-Stokes.
%
% Output arguments
% ----------------
%   TV, TVF : See notation.m

% Extract information from input
ndim = size(Qv1, 2)-1;
[nv1, ~, nq] = size(Qv1);
[~, ~, nqf, nf] = size(Qvf1);
[nv2, ~] = size(Qv2);

% Element basis
% Code me!
for i = 1:nv
    dof1 = (i - 1) * nv1 + 1: i * nv1;
    dof2 = (i - 1) * nv2 + 1: i * nv2;
    for k = 1:ndim + 1
        Tv(dof1, :, k, :) = pagemtimes(Qv1(i, k, :), [eye(nv1) zeros(nv1, nv2)]);
        Tvf(dof1, :, k, :, :) = pagemtimes(Qvf1(i, k, :, :), [eye(nv1) zeros(nv1, nv2)]);
        Tv(dof2, :, k, :) = pagemtimes(Qv2(i, k, :), [zeros(nv2, nv1) eye(nv2)]);
        Tvf(dof2, :, k, :, :) = pagemtimes(Qvf2(i, k, :, :), [zeros(nv2, nv1) eye(nv2)]);
    end
end

end