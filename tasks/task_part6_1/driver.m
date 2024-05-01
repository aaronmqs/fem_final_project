clear; clc; close all;

xlims = [0; 10]; ylims = [0; 1]; etype = 'hcube'; nel = [100; 10]; porder = 2; 
pltit = true;
[U, info] = solve_linelast_beam0(xlims, ylims, etype, nel, porder, pltit);