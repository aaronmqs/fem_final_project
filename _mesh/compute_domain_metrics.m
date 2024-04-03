function [v, c, sa] = compute_domain_metrics(transf_data, qrule)
%COMPUTE_DOMAIN_METRICS Compute the volume, centroid, and surface
%area of a domain (approximated) by the mesh described by TRANSF_DATA.
%
%Input arguments
%---------------
%  TRANSF_DATA, QRULE : See notation.m
%
%Output arguments
%----------------
%  V : number : Volume of domain
%
%  C : Array (NDIM,) : Centroid of domain
%
%  SA : number : Surface area of domain

% Extract relevant quantities
ndim = size(transf_data(1).xe, 1);
nf = size(transf_data(1).sigf, 2);
nelem = numel(transf_data);

% Initialize volume, centroid, surface area
c = zeros(ndim, 1);
v = 0.0; sa = 0.0;

% Code me!
[wq, wqf] = deal(qrule.wq, qrule.wqf);
v_fcn = ones(size(wq))';
s_fcn = ones(size(wqf))';

for e = 1:nelem
    % Volume
    transf_data_e = transf_data(e);
    [detG_e, sigf_e, e2bnd_e, xq_e] = deal(transf_data_e.detG, transf_data_e.sigf, transf_data_e.e2bnd, transf_data_e.xq);
    ve = v_fcn * (wq .* detG_e);
    v = v + ve;

    % Centroid
    c_fcn = xq_e;
    ce = c_fcn * (wq .* detG_e);
    c = c + ce;

    % Surface area
    sa_e = 0;
    for f = 1:nf
        if ~isnan(e2bnd_e(f))
            sa_e = sa_e + s_fcn * (wqf .* sigf_e(:, f));
        end
    end
    sa = sa + sa_e;
end

c = c / v;

end








