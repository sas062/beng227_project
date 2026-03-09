function [x,y,Rf] = circle_packing(Ri,r,thr,N)
    % r -> radius of cells (i.e smaller circle that is being packed into
    % larger circle)
    % Ri -> initial radius of islet (i.e larger circle that smaller circles are being
    % packed into)
    % thr -> number of failed placements required for algorithm to stop and
    % increase islet radius
    % N -> total number of cells desired

    % rng(1); % reproducible rand seed (can removed

    % store locations of smaller circles or 'cells'
    x = [];
    y = [];

    
    err = 0; % number of failed placements
    n = 0;   % number of circles packed

    Rcurr = Ri; % Current radius of islet

    while n < N
        if err >= thr
            Rcurr = Rcurr + r;
            err = 0;
        end
        r_cell = (Rcurr-r) * rand(1,1);
        phi = 2*pi * rand(1,1);
        x_cell = r_cell .* cos(phi);
        y_cell = r_cell .* sin(phi);

        b = check_overlap(x_cell,y_cell,x,y,r);
        if b == 1
            err = err+1;
            continue
        end
        err = 0;
        x = [x x_cell];
        y = [y y_cell];
        n = n+1;
    end
    Rf = Rcurr;
end

function b = check_overlap(xi,yi,x,y,r)
    b = 0;
    for i = 1:numel(x)
        dist = sqrt((xi-x(i))^2 + (yi-y(i))^2);
        if dist <= 2*r
            b = 1;
            break
        end
    end
end