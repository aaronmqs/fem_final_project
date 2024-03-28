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
if size(x, 1) ~= size(xk, 1); error("Incorrect dimension."); end
Q = eval_interp_simp_lagrange(xk, x);
[ndim, nx] = size(x);
tol = 1e-8;
for j = 1:nx
    if abs(sum(Q(:, 1, j)) - 1) > tol; error("The sum of the basis functions in a point must be 1."); end
    for k = 1:ndim
        if abs(sum(Q(:, 1 + k, j))) > tol; error("The sum of the derivatives of the basis functions in a point must be 0."); end
    end
end

% If this point is reached, the test was successful.
fprintf("Test passed.\n\n")

end