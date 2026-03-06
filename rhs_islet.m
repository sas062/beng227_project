function dydt = rhs_islet(t, y, p)

    % Indices
    ia = 1:p.Na;
    ib = p.Na + (1:p.Nb);
    iG = p.Na + p.Nb + (1:p.Nx*p.Ny);
    
    theta_a = y(ia);
    theta_b = y(ib);
    G = reshape(y(iG), p.Ny, p.Nx);
    
    insulin = insulin_secretion(theta_b, p);
    
    % Alpha-cell phase dynamics (Equation 1)
    dtheta_a = p.w_a + p.K_ba * insulin .* p.fra(theta_a);
    
    % Beta-cell phase dynamics (Equation 2)
    G_beta = interp2(p.X, p.Y, G, p.xb, p.yb, 'linear', 0);
    dtheta_b = p.w_b + p.K_ab * G_beta .* p.frb(theta_b);

    % Glucagon PDE dynamics (Equation 3)
    dGdt = glucagon_pde_rhs(t, G, theta_a, p);
    
    % Return column vector
    dydt = [dtheta_a; dtheta_b; dGdt(:)];

end
