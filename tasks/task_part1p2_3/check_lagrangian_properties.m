function check_lagrangian_properties(ndim, porder, x)
% This program checks the correctness of the implementation of the function
% that evaluates the Lagrangian basis functions and their derivatives for a
% d-dimensional simplex of order p.
%
% Test: properties of a Lagrangian basis (the sum of the basis functions in
% a point must be 1, while the sum of its derivatives must be zero)
%
% Example:
%
% check_lagrangian_properties(3, 5, [0.2 0.3; 0 0; 0.1 0.1])

xk = create_nodes_bndy_refdom_simp(ndim, porder);
nv = size(xk, 2);
if size(x, 1) ~= size(xk, 1); error("Incorrect dimension for evaluation points."); end
Q = eval_interp_simp_lagrange(xk, x);

% Checks
nx = size(x, 2);
tol = 1e-8;
Qk = eval_interp_simp_lagrange(xk, xk);
if norm(squeeze(Qk(:, 1, :)) - eye(nv), "inf") > tol; error("Lagrangian (nodal) property not satisfied."); end
if ~isempty(setdiff(find(abs(sum(squeeze(Q(:, 1, :)), 1) - 1) > tol), 1:nx))
    error("The sum of all basis functions should be equal to 1 for any given point.")
end
if ~isempty(setdiff(find((abs(sum(squeeze(Q(:, 2, :)), 1)) > tol)), 1:nx))
    error("The sum of the gradients of all basis functions should be equal to a zero vector for any given point.")
end

% If this point is reached, the test was successful.
fprintf("Test passed.\n\n")

end