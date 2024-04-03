% In this file the function that evaluates all relevant transformation
% quantities is tested.

clear; clc; close all;
tol = 1e-4;

%% Test for hcube
etype = 'hcube';
ndim = 2;
porder = 2;
nquad_per_dim = ceil((porder + 1) / 2) * 3; % TODO: Why so many qpoints necessary?

% Quadrature rules for ndim-dimensional hcube (\Omega_\square) and for
% (ndim-1)-dimensional hcube (\Gamma_\square)
qrule = create_qrule_gaussleg(etype, ndim, nquad_per_dim);
[wq, zq, wqf, rq] = deal(qrule.wq, qrule.zq, qrule.wqf, qrule.rq);

z = qrule.zq; r = qrule.rq;
lfcnsp = create_polysp_nodal(etype, ndim, porder, z, r);

xe = [0.0 0.6 1.2 0.2 0.8 1.6 0.0 0.6 1.2;
      0.0 0.4 0.2 0.8 1.0 1.0 1.6 1.4 1.8];
nv = (1 + porder)^ndim;
if size(xe, 2) ~= nv; error("Incorrect physical nodes."); end
Qv = lfcnsp.Qv;
Qvf = lfcnsp.Qvf;
r2z = lfcnsp.r2z;
f2v = lfcnsp.f2v;
[xq, detG, Gi, xqf, sigf, Gif] = eval_transf_quant_ndim(xe, Qv, Qvf, r2z, f2v);

v = ones(size(wq))';
V = v * (wq .* detG);
if abs(V - 1.6533) > tol; error("Incorrect volume."); end

c = xq;
C = (c * (wq .* detG)) / V;
if abs(C(1) - 0.8632) > tol || abs(C(2) - 0.9387) > tol; error("Incorrect centroid."); end

s = ones(size(wqf))';
S = 0;
for f = 1:size(sigf, 2)
    S = S + s * (wqf .* sigf(:, f));
end
if abs(S - 6.2791) > tol; error("Incorrect surface area."); end

% If this point is reached, the test was successful.
fprintf("Test for hypercube passed.\n\n")