function insulin = insulin_secretion(theta_b, p)

    W_tilde = exp( - ((p.xa - p.xb').^2 + (p.ya - p.yb').^2) / (2*p.lI^2) );
    W = W_tilde ./ sum(W_tilde, 2);
    
    insulin = W * p.fs(theta_b);

end
