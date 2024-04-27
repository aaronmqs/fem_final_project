clear; clc; close all;

cnt = 0; e = zeros(4, 5);
for porder = 1:5
    for k = 1:4
        k
        porder
        nel = 2^k; pltit = false; cnt = cnt + 1;
        [U, e(k, porder), info] = solve_pde0(nel, porder, pltit);
    end
    h = 1 ./ (porder * 2.^(1:4));
    ocr = log(e(1, porder) / e(end, porder)) / log(h(1) / h(end));
%     subplot(5, 1, porder)
    figure
    loglog(h, e(:, porder), 'Marker','square'); grid on;
    title("p = " + num2str(porder)+ ", Opt. conv. rate = " + num2str(ocr))
end