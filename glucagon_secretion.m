function S = glucagon_secretion(theta_a,p)

    S = zeros(p.Ny, p.Nx);
    
    for n = 1:p.Na
        amp = p.sG * p.fs(theta_a(n));
        phi_n = exp( - ( (p.X - p.xa(n)).^2 + (p.Y - p.ya(n)).^2 ) / (2*p.sigmaG^2) );
        phi_n = phi_n / (sum(phi_n(:)) * p.dx * p.dy); % normalize
        S = S + amp * phi_n;
    end
    
end