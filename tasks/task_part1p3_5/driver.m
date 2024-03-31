% This file checks if volumes and centroids for simplices and hypercubes in
% 1, 2 and 3 dimensions are evaluated correctly.

clear; clc; close all;

etype_ = ["simp" "hcube"]; porder = 1;
for etype = etype_    
    fprintf("Tests for element type = %s: \n\n", etype)
    for ndim = 1:3
        fprintf("Test for ndim = %d: \n", ndim);
        check_moments_refdom(etype, ndim, porder);
    end
end