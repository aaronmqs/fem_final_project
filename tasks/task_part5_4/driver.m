clc; close all; clear;

nel = [5 5 5]; porder = 2; pltit = true;
[U, info] = solve_linelptc_sclr_cube0(nel, porder, pltit);