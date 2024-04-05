function eval_interp_onedim_lagrange_pdeg(xk, porder)

p = porder;
pp1 = p + 1;
nx = 100; x = linspace(xk(1), xk(end), nx);
Q = eval_interp_onedim_lagrange(xk, x);

subplot(2,1,1)
for i = 1:pp1
    plot(x, squeeze(Q(i, 1, :)), LineWidth=2)
    hold on 
end
hold on
plot(xk, zeros(size(xk)), 'ob', LineWidth=1, MarkerSize=6)
grid on
ylabel("Basis functions")

subplot(2,1,2)
for i = 1:pp1
    plot(x, squeeze(Q(i, 2, :)), LineWidth=2)
    hold on 
end
plot(xk, zeros(size(xk)), 'ob', LineWidth=1, MarkerSize=6)
ylim([1.1*min(min(squeeze(Q(:, 2, :)))) 1.1*max(max(squeeze(Q(:, 2, :))))])
grid on
ylabel("Derivatives")
xlabel("x")

ttl = sprintf("%g-node Lagrangian functions (p = %g)", p+1, p);
sgtitle(ttl)

end