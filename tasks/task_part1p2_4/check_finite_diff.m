function check_finite_diff(ndim, porder, x)
% This program checks the correctness of the implementation of the function
% that evaluates the Lagrangian basis functions and their derivatives for a
% d-dimensional simplex of order p.
%
% Test: 
%
% Example:
%
% check_finite_diff(3, 5, [0.2 0.3; 0 0; 0.1 0.1])

xk = create_nodes_bndy_refdom_simp(ndim, porder); 
[ndim, nv] = size(xk);
if size(x, 1) ~= ndim; error("Incorrect dimension for evaluation points."); end
Q = eval_interp_simp_lagrange(xk, x);
e = eye(ndim); % canonical basis for R^ndim
eps = 1e-7;
nx = size(x, 2);

% Finite difference approximation
% Fd(i, j, k) is the finite difference approximation to the derivative of 
% the ith basis function with respect to the jth coordinate evaluated at 
% the kth point. We want to compare this with Q(i, 1 + j, k).
Fd = zeros(size(Q(:, 2:end, :)));
for j = 1:ndim
    Ej = repmat(e(:, j), [1, nx]);
    Qp = eval_interp_simp_lagrange(xk, x + eps * Ej);
    Qm = eval_interp_simp_lagrange(xk, x - eps * Ej);
    Fd(:, j, :) = (Qp(:, 1, :) - Qm(:, 1, :)) / (2 * eps);
end

% Finite difference test
tol = 10 * eps;
for i = nv
    for k = 1:nx
        if norm(Q(i, 2:end, k) - Fd(i, :, k)) > tol
            error("The finite difference test failed.")
        end
    end
end

% If this point is reached, the test was successful.
fprintf("Test passed.\n\n")

end







