function dydt = rhs_islet(t, y, p)

    % Indices
    ia = 1:p.Na;
    ib = p.Na + (1:p.Nb);
    iG = p.Na + p.Nb + (1:p.Nx*p.Ny);
    %iI = p.Na + p.Nb + p.Nx*p.Ny + (1:p.Nx*p.Ny);
    
    theta_a = y(ia);
    theta_b = y(ib);
    G = reshape(y(iG), p.Ny, p.Nx);

    % Boundary Conditions -> Perfectly absorbing
    [n,m] = size(G);
    G(1:n,1) = 0;
    G(1:n,m) = 0;
    G(1,1:m) = 0;
    G(n,1:m) = 0;

    % Alpha-cell phase dynamics (Equation 1)
    insulin = insulin_secretion(theta_b, p);
    dtheta_a = p.w_a + p.K_ba * insulin .* p.fra(theta_a);
    
    % Beta-cell phase dynamics (Equation 2)
    G_beta = interp2(p.X, p.Y, G, p.xb, p.yb, 'linear', 0);
    dtheta_b = p.w_b + p.K_ab * p.fGLP1R(G_beta) .* p.frb(theta_b);

    % Beta-beta interactions (cx-36 coupling)
    dtheta_bb = zeros(1,p.Nb);
    for i = 1:p.Nb
        NNi = p.NN(i,:);
        [~,cells,~] = find(NNi);
        if numel(cells) == 0
            continue
        end
        dtheta_bb_i = 0;
        for j = 1:numel(cells)
            dtheta_bb_i = dtheta_bb_i + p.K_bb*sin(theta_b(j)-theta_b(i)); % Kuramato model
        end    
        dtheta_bb(i) = dtheta_bb_i/numel(cells);
    end
    dtheta_b = dtheta_b + dtheta_bb';

    % Glucagon PDE dynamics (Equation 3)
    dGdt = glucagon_pde_rhs(t, G, theta_a, p);
    
    % Return column vector
    dydt = [dtheta_a; dtheta_b; dGdt(:)];
    disp(t)
end
