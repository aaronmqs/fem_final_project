clc; clear; close all;

nref = 0; porder = 2; pltit = true;
[U, info] = solve_linelptc_sclr_batman0(nref, porder, pltit);

save("U", "U")
save("info", "info")