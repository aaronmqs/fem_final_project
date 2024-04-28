function [] = fdtest(f, x)

[~, dF] = f(x);
d = length(x);
dF_fd = zeros(size(dF));
eps = 1e-7; e = eye(d);
for j = 1:d
    Fp = f(x + eps * e(:, j));
    Fm = f(x - eps * e(:, j));
    dF_fd(:, j) = (Fp - Fm) / (2 * eps);
end
if norm(dF_fd - dF) > 1e-6 ; error("Incorrect jacobian function."); end

fprintf("Test passed");

end