clc; clear; close all;

etype = 'hcube'; nel = [20 20]; porder = 2; pltit = true;
[U, E, info] = solve_linelptc_sclr_disk0(etype, nel, porder, pltit);