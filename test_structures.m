close all; clear all; clc;

% This script tests the circle packing algorithm

N = 100;   % number of cells
Ri = 10;   % um
thr = 100; % number of errors before increasing islet radius
r = 5;     % radius of cell
rng(1);


[x,y,Rf] = circle_packing(Ri,r,thr,N);

% Plot islet structure
figure(1);
axis equal tight
for i = 1:N
    rectangle('Position',[x(i)-r y(i)-r 2*r 2*r],'Curvature',[1,1],'EdgeColor','k');
end
rectangle('Position',[-Rf -Rf 2*Rf 2*Rf],'Curvature',[1,1],'EdgeColor','k')


%% Now lets get connectivity matrix

Na = floor(N*0.15);
Nb = N-Na;
arrangement = "toward_outer";
%arrangement = "mixed";
%arrangement = "outer";
[xa,ya,xb,yb] = assign_cell_type(x,y,Na,Nb,arrangement);

thr = 20;
NN = get_NN(xb,yb,thr);

figure(4)
for i=1:Na
    rectangle('Position',[xa(i)-r ya(i)-r 2*r 2*r],'Curvature',[1,1],'FaceColor','g','EdgeColor','g');
end

for i=1:Nb
    rectangle('Position',[xb(i)-r yb(i)-r 2*r 2*r],'Curvature',[1,1],'FaceColor', 'r', 'EdgeColor','r');
end
rectangle('Position',[-Rf -Rf 2*Rf 2*Rf],'Curvature',[1,1],'EdgeColor','k')

for i = 1:Nb
    xi = xb(i);
    yi = yb(i);
    NNi = NN(i,:);
    [~,cells,~] = find(NNi);
    for j = 1:numel(cells)
        xj = xb(cells(j));
        yj = yb(cells(j));
        line([xi xj], [yi yj],'Color', 'black','LineWidth',1)
    end    
end

structural_connections = sum(NN);