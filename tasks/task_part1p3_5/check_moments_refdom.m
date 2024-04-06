function lfcnsp = check_moments_refdom(etype, ndim, porder)
% This function checks the quadrature formulas provided by computing the
% moments corresponding to Volume and Centroid of a master domain (etype = simplex or hypercube)

% A 1-dimensional polynomial of degree p = 2 * nq - 1 can be integrated exactly with nq
% quadrature points by using a Gaussian quadrature rule. This can be
% extended to the d-dimensional by a tensor product of 1d quadrature rules,
% even if it's not optimal (a smaller number of points could possibly be
% used).
if strcmp(etype, "hcube")
    nquad_per_dim = 1;
elseif strcmp(etype, "simp")
    % It was necessary to use one more quadrature point for simplices.
    % Check case: check_moments_refdom("simp", 2, 1)
    nquad_per_dim = 2;
end

% Quadrature rules for ndim-dimensional simplex (\Omega_\square) and for
% (ndim-1)-dimensional simplex (\Gamma_\square)
qrule = create_qrule_gaussleg(etype, ndim, nquad_per_dim);
[wq, zq, wqf, rq] = deal(qrule.wq, qrule.zq, qrule.wqf, qrule.rq);

% Volume of ndim-dimensional simplex V(\Omega_\square)
vo = ones(size(wq))';
Vo = vo * wq;

% Centroid of ndim-dimensional simplex V(\Omega_\square)
co = zq;
Co = (co * wq) / Vo;

% Volume of (ndim-1)-dimensional simplex V(\Gamma_\square)
vg = ones(size(wqf))';
Vg = vg * wqf;

% Centroid of ndim-dimensional simplex V(\Omega_\square)
cg = rq;
Cg = (cg * wqf) / Vg;

% Local function space
z = qrule.zq; r = qrule.rq;
if nargout > 0; lfcnsp = create_polysp_nodal(etype, ndim, porder, z, r); end

tol = 1e-8;
% Checks volumes
if strcmp(etype, "simp")
    if abs(Vo - 1 / factorial(ndim)) > tol; error("Incorrect moment."); end
    if abs(Vg - 1 / factorial(ndim - 1)) > tol; error("Incorrect moment."); end
elseif strcmp(etype, "hcube")
    if abs(Vo - 2^ndim) > tol; error("Incorrect moment."); end
    if abs(Vg - 2^(ndim - 1)) > tol; error("Incorrect moment."); end
else
    warning("Case not considered in the tests.")
end

% Checks centroids
for k = 1:ndim
    if strcmp(etype, "simp")
        if abs(Co(k) - 1 / (ndim + 1)) > tol; error("Incorrect moment."); end
        if k < ndim && abs(Cg(k) - 1 / ndim) > tol; error("Incorrect moment."); end
    elseif strcmp(etype, "hcube")
        if abs(Co(k)) > tol; error("Incorrect moment."); end
        if k < ndim && abs(Cg(k)) > tol; error("Incorrect moment."); end
    else
        warning("Case not considered in the tests.")
    end    
end

% If this point is reached, the test was successful.
fprintf("Test passed.\n\n")

end




