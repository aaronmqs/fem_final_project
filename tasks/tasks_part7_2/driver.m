clc; clear; close all;

rho = 1; L = 1; etype = 'simp'; porder = 2; 
pltit = true; U = 1; U0 = [];

for Re = 100:100:1200
    nu = U * L / Re;
    U0 = solve_ins_nd0(rho, nu, 0, porder, false, U0);
end

Re = 1300;
nu = U * L / Re;
U0 = solve_ins_nd0(rho, nu, 0, porder, true, U0);