function dGdt = glucagon_pde_rhs(t, G, theta_a, p)

    % Discrete Laplacian
    L = del2(G, p.dx, p.dy);
    
    S = glucagon_secretion(theta_a,p);
    
    % PDE right-hand side
    dGdt = p.DG * L + S - p.tau * G;

end
