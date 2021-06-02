function [h1, h2] = plot_dir3 (vX, vY, vZ)
%function [h1, h2] = plot_dir3 (vX, vY, vZ)
%Plotting x-y variables with direction indicating vector to the next element.
%Example
%   vX = linspace(0,2*pi, 10)';
%   vY = sin (vX);
%   vZ = cos (vX);
%   plot_dir3(vX, vY, vZ);
colorR = [0.86 0.30 0.43];
rMag = 0.5;

% Length of vector
lenTime = length(vX);

% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;

% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);
% X coordinates of tails of arrows
vZQ0 = vZ(vSelect0, 1);

% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);
% Z coordinates of heads of arrows
vZQ1 = vZ(vSelect1, 1);

% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;
vPz = (vZQ1 - vZQ0) * rMag;

% make plot 
h1 = plot3 (vX, vY, vZ, 'c-.','LineWidth',0.5,'MarkerSize',3); hold on;
set(h1,'Marker','.');
set(h1,'MarkerEdgeColor',colorR);
set(h1,'MarkerFaceColor','y');
% add arrows 
% h2 = quiver3 (vXQ0, vYQ0, vZQ0, vPx, vPy, vPz, 0, '-.c'); grid on; hold off
axis equal