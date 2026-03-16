function S = glucagon_secretion(theta_a,p)

    S = zeros(p.Ny, p.Nx);
    
    for n = 1:p.Na
        V = 250*(1e-15); % L x W x H  = 5 x 5 x 10 = 250 um, conversion um^3 -> L
        amp = p.sG_G * p.fs(theta_a(n));
        phi_n = exp( - ( (p.X - p.xa(n)).^2 + (p.Y - p.ya(n)).^2 ) / (2*p.sigmaG^2) );
        phi_n = phi_n / (sum(phi_n(:)) * p.dx * p.dy); % normalize
        S = S + amp * phi_n/V*1e-18*1e9;% Conversion amol to mol and M -> nM
    end
    
end