clear; clc; close all;

c = [0; 0]; r1 = 1; r2 = 2; h = 10; nel = [20; 4; 15]; porder = 2;
pltit = true;
[U, info] = solve_linelast_hcyl0(c, r1, r2, h, nel, porder, pltit);