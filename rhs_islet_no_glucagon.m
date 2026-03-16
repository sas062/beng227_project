function dydt = rhs_islet_no_glucagon(t, y, p)

    % Indices
    ia = 1:p.Na;
    ib = p.Na + (1:p.Nb);
    
    theta_a = y(ia);
    theta_b = y(ib);
    
    %insulin = insulin_secretion(theta_b, p);
    
    % Alpha-cell phase dynamics (Equation 1)
    dtheta_a = ones(p.Na,1)*p.w_a; %+ p.K_ba * insulin .* p.fra(theta_a);
    
    % Beta-cell phase dynamics (Equation 2)
    dtheta_b = p.w_b; 

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
            dtheta_bb_i = dtheta_bb_i + p.K_bb*sin(theta_b(j)-theta_b(i)); % kuramato model
        end    
        dtheta_bb(i) = dtheta_bb_i/numel(cells);
    end
    dtheta_b = dtheta_b + dtheta_bb';


    
    % Return column vector
    dydt = [dtheta_a; dtheta_b];

end
