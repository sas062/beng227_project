function [xa,ya,xb,yb] = assign_cell_type(x,y,Na,Nb,arrangement)
    % This code assigns cell type to the islet structure (i.e. determines
    % which cells are alpha cells and which are beta cells

    % [x,y] -> locations of cells
    % Na -> number of alpha-cells
    % Nb -> number of beta-cells
    % arrangement -> mixed, outer, or toward_outer
    
    %rng(1); % reproducible rand seed
    N = Na+Nb; % total number of cells

    if arrangement == "mixed" % will mix alpha and beta cells together
        indices_alpha = randperm(N,Na);
        xa = x(indices_alpha);
        ya = y(indices_alpha);
        
        %indices_beta = 1:N;
        %indices_beta(indices_alpha) = [];
        xb = x(setdiff(1:N,indices_alpha));
        yb = y(setdiff(1:N,indices_alpha));
    end

    if arrangement == "outer" % will make alpha cells located on the periphery
        for i = 1:N
            dist(i) = sqrt(x(i)^2 + y(i)^2);
        end
        [~,I] = sort(dist,'descend');
        xa = x(I(1:Na));
        xb = x(I((Na+1):end));

        ya = y(I(1:Na));
        yb = y(I((Na+1):end));
    end

    if arrangement == "toward_outer" % This will arrange alpha-cells towards the periphery, but not completely at the periphery
        for i = 1:N
            dist(i) = sqrt(x(i)^2 + y(i)^2);
        end
        [~,I] = sort(dist,'descend');
        potential_alpha = I(1:floor(N/2));
        indices_alpha = potential_alpha(randperm(floor(N/2),Na));
        disp(indices_alpha)

        xa = x(indices_alpha);
        ya = y(indices_alpha);

        % indices_beta = [];
        % for i = 1:N
        %     if numel(find(indices_alpha == i)) == 0
        %         indices_beta = [indices_beta i];
        %     end
        % indices_beta = I;
        % indices_beta(indices_alpha) = [];
        % 
        % xb = x(indices_beta);
        % yb = y(indices_beta);
        
        xb = x(setdiff(1:N,indices_alpha));
        yb = y(setdiff(1:N,indices_alpha));

    end

end