function [Tv, Tvf] = create_elem_basis(nvar, Qv, Qvf)
%CREATE_ELEM_BASIS Create element basis evaluated at point throughout
%volume (TV) and on each face (TVF) from basis of local function space
%evaluated at corresponding points (QV, QVF).
%
% Input arguments
% ---------------
%   NVAR, QV, QVF : See notation.m
%
% Output arguments
% ----------------
%   TV, TVF : See notation.m

% Extract information from input
ndim = size(Qv, 2)-1;
[nv, ~, nq] = size(Qv);
[~, ~, nqf, nf] = size(Qvf);

% Element basis
% Code me!

% Preallocate
ndof_per_elem = nvar * nv;
Tv = zeros(ndof_per_elem, nvar, ndim+1, nq);
Tvf = zeros(ndof_per_elem, nvar, ndim+1, nqf, nf);

% Evaluate Tv
for q = 1:nq
    for i = 1:nv
        dof = (i - 1) * nvar + 1: i * nvar;
        for k = 1:ndim+1
            Tv(dof, :, k, q) = Qv(i, k, q) * eye(nvar);
        end
    end
end

% Evaluate Tvf
for q = 1:nqf
    for i = 1:nv
        dof = (i - 1) * nvar + 1: i * nvar;
        for k = 1:ndim+1
            for f = 1:nf
                Tvf(dof, :, k, q, f) = Qvf(i, k, q, f) * eye(nvar);
            end
        end
    end
end

end







