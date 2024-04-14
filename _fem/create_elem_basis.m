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
for i = 1:nv
    dof = (i - 1) * nvar : i * (nvar + 1) - 1;
    Tvf(dof, :, :, :, :) = kron(Qvf(i, :, :, :), eye(nvar));
    Tv(dof, :, :, :) = kron(Qv(i, :, :), eye(nvar));
end

% for i = 1:nv
%     dof = (i - 1) * nvar : i * (nvar + 1) - 1;
%     for k = 1:ndim + 1
%         for l = 1:nq
%             Tv(dof, :, k, l) = Qv(i, k, l) * eye(nvar);
%         end
%         for l = 1:nqf
%             for f = 1:nf
%                 Tvf(dof, :, k, l, f) = Qvf(i, k, l, f) * eye(nvar);
%             end
%         end
%     end
% end

end











