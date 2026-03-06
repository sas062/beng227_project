function dGdt = glucagon_pde_rhs(t, G, theta_a, p)

    % Discrete Laplacian
    L = del2(G, p.dx, p.dy);

    % Explicit 5-pt stencil implementation
    % L = zeros(size(G));
    % L(2:end-1,2:end-1) = ...
    % (G(2:end-1,3:end) - 2*G(2:end-1,2:end-1) + G(2:end-1,1:end-2)) / p.dx^2 + ...
    % (G(3:end,2:end-1) - 2*G(2:end-1,2:end-1) + G(1:end-2,2:end-1)) / p.dy^2;

    S = glucagon_secretion(theta_a,p);
    
    % PDE right-hand side
    dGdt = p.DG * L + S - p.tau * G;

end
