clear; clc; close all;

kmax = 4;
pmax = 5;
ocr = zeros(pmax, 1);
e = zeros(4, 5);
for porder = 1:5
    for k = 1:kmax
        nel = 2^k; pltit = false;
        [U, e(k, porder), info] = solve_pde0(nel, porder, pltit);
    end
    h = 1 ./ (2.^(1:kmax));
    if porder < 5
        ocr(porder) = log(e(1, porder) / e(end, porder)) / log(h(1) / h(end));
    else
        ocr(porder) = log(e(1, porder) / e(end-1, porder)) / log(h(1) / h(end-1));
    end
%     subplot(5, 1, porder)
    loglog(h, e(:, porder), 'Marker','square'); grid on; hold on;
%     title("p = " + num2str(porder)+ ", Opt. conv. rate = " + num2str(ocr(porder)))
    fprintf("p = " + num2str(porder)+ ", Opt. conv. rate = " + num2str(ocr(porder)) + "\n\n")
end

legend("p1", "p2", "p3", "p4", "p5")