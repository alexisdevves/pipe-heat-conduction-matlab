function plot_mesh(rI, rO, Nr, Nth)
% PLOT_MESH  Draw the polar grid of the pipe wall and number the nodes.
%
%   rI, rO : inner / outer radius [m]
%   Nr     : number of radial nodes (rings)
%   Nth    : number of angular nodes (spokes) on the half section 0..pi
%
%   The full circle is drawn, but nodes are numbered only on the solved half
%   (left half, starting from the top). Node id = (g-1)*Nr + ir.

nIntervals   = Nr - 1;                    % number of radial intervals
dr           = (rO - rI)/nIntervals;      % radial step
angleStepDeg = 180/(Nth - 1);             % angular step [deg]

t = linspace(0, 2*pi, 400);

light_figure;
hold on; grid on;

%% ---- Fill the first radial layer (between rI and rI+dr) ----
x_in  = rI*cos(t);        y_in  = rI*sin(t);
x_out = (rI + dr)*cos(t); y_out = (rI + dr)*sin(t);
fill([x_in fliplr(x_out)], [y_in fliplr(y_out)], 'r', ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');

%% ---- Rings ----
for k = 0:nIntervals
    rk = rI + k*dr;
    plot(rk*cos(t), rk*sin(t), 'k', 'LineWidth', 2);
end

%% ---- Spokes over the full circle ----
angles_deg = (0:2*(Nth - 1) - 1)*angleStepDeg;
for ang = angles_deg
    a = deg2rad(ang);
    plot([rI*cos(a) rO*cos(a)], [rI*sin(a) rO*sin(a)], 'k', 'LineWidth', 1.5);
end

%% ---- Nodes and numbering (half section: 90..270 deg, top to bottom) ----
node_angles_deg = 90 + (0:Nth-1)*angleStepDeg;
offset = 0.1*(rO - rI);       % shift the label away from the grid lines
nodeID = 1;

for ang = node_angles_deg
    a = deg2rad(ang);

    for k = 0:nIntervals
        rk = rI + k*dr;
        x = rk*cos(a);
        y = rk*sin(a);

        plot(x, y, 'k.', 'MarkerSize', 14);

        text(x - offset*sin(a), y + offset*cos(a), num2str(nodeID), ...
            'Color', 'k', 'FontSize', 9, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

        nodeID = nodeID + 1;
    end
end

axis equal;
axis off;
xlim(1.12*[-rO rO]);
ylim(1.12*[-rO rO]);
hold off;
end
