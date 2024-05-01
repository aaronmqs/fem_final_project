clc; close all; clear;

nel = [10 10 10]; porder = 2; pltit = true;
[U, info] = solve_linelptc_sclr_cube0(nel, porder, pltit);