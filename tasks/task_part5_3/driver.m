clear; clc; close all;

nref = 0; porder = 2; pltit = true;
[U, info] = solve_linelptc_sclr_nd0(nref, porder, pltit);