function insulin = insulin_secretion(theta_b, p)

    W_tilde = exp( - ((p.xa - p.xb').^2 + (p.ya - p.yb').^2) / (2*p.lI^2) );
    W = W_tilde ./ sum(W_tilde, 2);

    insulin = W * p.fs(theta_b);

    % S = zeros(p.Ny, p.Nx);
    % 
    % for n = 1:p.Nb
    %     amp = p.sG_I * p.fs(theta_b(n));
    %     phi_n = exp( - ( (p.X - p.xb(n)).^2 + (p.Y - p.yb(n)).^2 ) / (2*p.sigmaG^2) );
    %     phi_n = phi_n / (sum(phi_n(:)) * p.dx * p.dy); % normalize
    %     S = S + amp * phi_n;
    % end
    % 
    % insulin = S;

end
