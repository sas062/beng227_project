function NN = get_NN(xb,yb,thr)
    % This code creates the connectivity matrix for beta-beta interactions
    % (i.e. cx-36 gap junction coupling)
    % [xb, yb] -> locations of beta-cells
    % thr -> distance threshold (um) for cells to be considered coupled
    Nb = numel(xb);

    NN = zeros(Nb);

    for i = 1:Nb
        for j = 1:Nb
            if i ~= j
                dist_ij = sqrt((xb(i)-xb(j))^2+(yb(i)-yb(j))^2);
                if dist_ij < thr
                    NN(i,j) = 1;
                    NN(j,i) = 1;
                end
            end
    
        end
    end

end