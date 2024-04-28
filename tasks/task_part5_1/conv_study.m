clear; clc; close all;

pmax = 3;
kmax = 4;
pltit = false;

% Simplex
etype = 'simp';
Es = zeros(kmax, pmax); h = Es; ocr = zeros(pmax, 1);
for porder = 1:pmax
    for k = 1:kmax
        nel = 2^k * [2 2];

        % Compute error
        [~, Es(k, porder)] = solve_linelptc_sclr_disk0(etype, nel, porder, pltit);

        % Compute element size
        msh = create_mesh_hsphere(etype, [0; 0], 1, nel, porder);
        v = measure_domain(msh, porder, msh.ndim, etype);
        h(k, porder) = sqrt(v / size(msh.e2vcg, 2));
    end
    ocr(porder) = log(Es(1, porder) / Es(end, porder)) / log(h(1) / h(end));
    loglog(h, Es(:, porder), 'Marker','square'); grid on; hold on;
    fprintf("p = " + num2str(porder)+ ", Opt. conv. rate = " + num2str(ocr(porder)) + "\n\n")
end

legend("p1", "p2", "p3")

save("Es", "Es")
save("h", "h")