function R = network_analysis(Kavg_desired,theta_b,fspan)
    st = fspan(1);
    ed = fspan(2);
    sc = 0.01; % Stopping criterion of 1 %;
    [frames,numCells] = size(theta_b);

    rho = corr(theta_b(st:ed,:));

    % These are your starting bounds for R. The desired R must be
    % within the range [R_low R_high].This means you must choose a value 
    R_low = 0.00; 
    R_high = 1;
    
    R = (R_low+R_high)/2; % Starting guess

    Kavg_low = mean(calc_links(numCells,R_low,rho));
    Kavg_high = mean(calc_links(numCells,R_high,rho));

    Kavg = mean(calc_links(numCells,R,rho));

% This basically uses the bisection root finding method to reach a desired
% Kavg

    while abs(Kavg - Kavg_desired)/Kavg_desired*100 > sc % Calculate the percent difference to desired kavg and compare to stopping criterion
        disp(R);
        diff_Kavg = Kavg - Kavg_desired;
        if diff_Kavg > 0
            R_low = R;
        else
            R_high = R;
        end
    
        R = (R_low+R_high)/2;
        Kavg_low = mean(calc_links(numCells,R_low,rho));
        Kavg_high = mean(calc_links(numCells,R_high,rho));
    
        Kavg = mean(calc_links(numCells,R,rho));
    end


    % Plot Degree Distribution
    
    numLinks = calc_links(numCells,R,rho);
    numLinksMax = max(numLinks);
    edges = 0:(numLinksMax+1);
    figure;
    hold on;
    title("Number of Cells vs. Number of Links");
    xlabel("Number of Links");
    ylabel("Number of Cells");
    histogram(numLinks,edges)
    
end

function numLinks = calc_links(n, threshold, rho)
    % NUMLINKS => 1 by n vector where each column,i, represents cell i's
    % number of links
    % INPUT => your calcium trace in whic columns represent the cell number
    % and rows represent the timepoints
    % THRESHOLD => the threshold for defining a link
    % This uses the builtin corr() function to calcualte the number of
    % links, which is equivalent to Dr. Kravet Network Analysis code
numLinks = zeros(1,n);
for i = 1:n
    for j = 1:n
        if i == j
            continue
        end
        if rho(i,j) > threshold
            numLinks(i) = numLinks(i)+1;
        end
    end
end



end