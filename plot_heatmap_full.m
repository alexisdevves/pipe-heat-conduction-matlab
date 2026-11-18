function plot_heatmap_full(T, Nr, Nth, rI, rO, theta, iM)
% PLOT_HEATMAP_FULL  Full 360-degree polar heatmap with the grid and node numbers.
%
%   T      : solution vector (Nr*Nth x 1)
%   Nr,Nth : number of radial / angular nodes
%   rI,rO  : inner / outer radius [m]
%   theta  : spoke angles on the half section (1 x Nth, 0..pi)
%   iM     : index of the ring on the material interface (drawn thicker)

[X, Y, Tfull, r_nodes] = interpolate_full_field(T, Nr, Nth, rI, rO, theta, iM);

% Spoke angles on the plot: theta = 0 at the top, solved half on the left (90..270 deg)
theta_half_plot = theta + pi/2;

%% Heatmap
light_figure;
pcolor(X, Y, Tfull);
shading interp;
axis equal;
axis off;
colorbar;
title('Temperature field');
hold on;

%% Grid: rings (interface ring thicker)
tt = linspace(0, 2*pi, 800);
for ir = 1:Nr
    rr = r_nodes(ir);
    if ir == iM
        lw = 2.0;
    else
        lw = 1.0;
    end
    plot(rr*cos(tt), rr*sin(tt), 'k-', 'LineWidth', lw);
end

%% Grid: spokes (solved half and its mirror image)
for g = 1:Nth
    th    = theta_half_plot(g);
    th_op = th + pi;

    plot([rI*cos(th)    rO*cos(th)],    [rI*sin(th)    rO*sin(th)],    'k-', 'LineWidth', 1.3);
    plot([rI*cos(th_op) rO*cos(th_op)], [rI*sin(th_op) rO*sin(th_op)], 'k-', 'LineWidth', 1.3);
end

%% Nodes and numbering (solved half only)
nodeID = 1;
offset = 0.1*(rO - rI);

for g = 1:Nth
    th = theta_half_plot(g);

    for ir = 1:Nr
        rr = r_nodes(ir);
        x = rr*cos(th);
        y = rr*sin(th);

        plot(x, y, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');

        text(x - offset*sin(th), y + offset*cos(th), num2str(nodeID), ...
            'Color', 'k', 'FontSize', 8, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

        nodeID = nodeID + 1;
    end
end

%% Wall boundaries
plot(rI*cos(tt), rI*sin(tt), 'k', 'LineWidth', 2.2);
plot(rO*cos(tt), rO*sin(tt), 'k', 'LineWidth', 2.2);

% Margin so the labels of the outer ring are not clipped
xlim(1.1*[-rO rO]);
ylim(1.1*[-rO rO]);
hold off;
end
