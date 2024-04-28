clear; clc; close all;

pmax = 3;
kmax = 4;
pltit = false;
load("Eh.mat"); load("Es.mat"); load("hh.mat"); load("hs.mat");

%% Simplex
ocr = zeros(pmax, 1);
for porder = 1:pmax
    ocr(porder) = log(Es(1, porder) / Es(end, porder)) / log(hs(1, porder) / hs(end, porder));
    loglog(hs(:, porder), Es(:, porder), 'Marker','square'); grid on; hold on;
    fprintf("\n p = " + num2str(porder)+ ", Opt. conv. rate = " + num2str(ocr(porder)) + "\n\n")
end

legend("p1", "p2", "p3")

%% Hypercube
figure
ocr = zeros(pmax, 1);
for porder = 1:pmax
    ocr(porder) = log(Eh(1, porder) / Eh(end, porder)) / log(hh(1, porder) / hh(end, porder));
    loglog(hh(:, porder), Eh(:, porder), 'Marker','square'); grid on; hold on;
    fprintf("p = " + num2str(porder)+ ", Opt. conv. rate = " + num2str(ocr(porder)) + "\n\n")
end

legend("p1", "p2", "p3")