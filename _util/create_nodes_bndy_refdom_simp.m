function [zk, f2v, N] = create_nodes_bndy_refdom_simp(ndim, porder)
%CREATE_NODES_BNDY_REFDOM_SIMP Create nodal distribution and boundary of
%NDIM-dimensional simplex element of order PORDER.
%
%Input arguments
%---------------
%   NDIM, PORDER : See notation.m
%
%Output arguments
%----------------
%   ZK, F2V, N : See notation.m

% Treat 0-dimensional case (boundary of 1D element) as special case
if ndim == 0
    zk = zeros(0, 1); f2v = []; N = [];
    return;
end

% Extract information from input
nf = ndim+1;

% Create nodal distribution for simplex
zk = tensprod_vector_from_onedim_unif(linspace(0, 1, porder+1), ndim);
zk = zk(:, sum(zk, 1)<=1.0);
if porder == 0, zk = [0.5; 0.5]; end

% Create a helper matrix as a tensor product of 1:porder+1 that will assist
% in extract face information for the faces of the unit simplex that
% parallel faces of the regular hcube
shp = (porder+1)*ones(1, ndim);
M = nan(shp);
cnt = 0;
nv_dm1 = 0;
for i=1:numel(M)
    idx = cell(1, ndim);
    [idx{:}] = ind2sub(shp, i);
    idx = cell2mat(idx);
    if sum(idx-1)>porder, continue; end
    if sum(idx-1)==porder, nv_dm1=nv_dm1+1; end
    cnt = cnt + 1;
    M(i) = cnt;
end

% Create mapping from face to nodes of element, first ndim faces (aligned
% with coordinate axes, i.e., parallel to faces of regular hcube)
f2v = zeros(nv_dm1, nf);
for i = 1:ndim
    idx = cell(1, ndim);
    [idx{:}] = deal(1:porder+1);
    
    idx{i} = 1;
    nodes = sort(squeeze(M(idx{:})));
    f2v(:, i) = nodes(~isnan(nodes));
end

% Create mapping from face to nodes for the inclined face
cnt=0;
for i=1:numel(M)
    idx = cell(1, ndim);
    [idx{:}] = ind2sub(shp, i);
    idx = cell2mat(idx);
    if sum(idx-1)==porder, cnt=cnt+1; f2v(cnt, end) = M(i); end
end

% Create unit normals for each face
N = zeros(ndim, nf);
for i = 1:ndim
    N(i, i) = -1;
end
N(:, end) = 1/sqrt(ndim);

end