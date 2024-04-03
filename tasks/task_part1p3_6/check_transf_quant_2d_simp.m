% In this file the function that evaluates all relevant transformation
% quantities is tested.

clear; clc; close all;
tol = 1e-4;

%% Test for simplex
etype = 'simp';
ndim = 2;
porder = 2;
nquad_per_dim = ceil((porder + 1) / 2) + 1; % More points necessary than "expected" because the integrand "suffers" two different transformations: one from the idealized hypercube to the idealized simplex and another one from the idealized simplex to the physical one.

% Quadrature rules for ndim-dimensional simplex (\Omega_\square) and for
% (ndim-1)-dimensional simplex (\Gamma_\square)
qrule = create_qrule_gaussleg(etype, ndim, nquad_per_dim);
[wq, zq, wqf, rq] = deal(qrule.wq, qrule.zq, qrule.wqf, qrule.rq);

z = qrule.zq; r = qrule.rq;
lfcnsp = create_polysp_nodal(etype, ndim, porder, z, r);

xe = [0 0.5 1.0 -0.3 0.3 -0.25;
      0 0.15 0.7 0.5 0.75 1.2];
nv = nchoosek(ndim + porder, ndim);
if size(xe, 2) ~= nv; error("Incorrect physical nodes."); end
Qv = lfcnsp.Qv;
Qvf = lfcnsp.Qvf;
r2z = lfcnsp.r2z;
f2v = lfcnsp.f2v;
[xq, detG, Gi, xqf, sigf, Gif] = eval_transf_quant_ndim(xe, Qv, Qvf, r2z, f2v);

v = ones(size(wq))';
V = v * (wq .* detG);
if abs(V - 0.7858) > tol; error("Incorrect volume."); end

c = xq;
C = (c * (wq .* detG)) / V;
if abs(C(1) - 0.1805) > tol || abs(C(2) - 0.4993) > tol; error("Incorrect centroid."); end

s = ones(size(wqf))';
S = 0;
for f = 1:size(sigf, 2)
    S = S + s * (wqf .* sigf(:, f));
end
if abs(S - 4.0158) > tol; error("Incorrect surface area."); end

% If this point is reached, the test was successful.
fprintf("Test for simplex passed.\n\n")