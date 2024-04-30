clear; clc; close all;

rho = 1; L = 1; etype = 'simp'; nel = [20; 20]; porder = 2; 
pltit = true; U = 1; U0 = [];

for Re = 100:100:1900
    nu = U * L / Re;
    U0 = solve_ins_ldc0(rho, nu, L, etype, nel, porder, false, U0);    
end

Re = 2000; nu = U * L / Re;
[U0, info] = solve_ins_ldc0(rho, nu, L, etype, nel, porder, pltit, U0);