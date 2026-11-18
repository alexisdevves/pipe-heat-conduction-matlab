function plot_node_temperatures(T, Nr, Nth, rI, rO, theta, rings)
% PLOT_NODE_TEMPERATURES  Write the nodal temperatures at the node positions,
% only for the selected rings.
%
%   T      : solution vector (Nr*Nth x 1), node id = (g-1)*Nr + ir
%   Nr,Nth : number of radial / angular nodes
%   rI,rO  : inner / outer radius [m]
%   theta  : spoke angles on the half section (1 x Nth, 0..pi)
%   rings  : ring indices to show, e.g. [1 2 4 6]

if any(rings < 1 | rings > Nr)
    error('rings must contain indices between 1 and Nr = %d.', Nr);
end

dr = (rO - rI)/(Nr - 1);
r  = rI + (0:Nr-1)*dr;              % ring radii

% Plot angles: theta = 0 at the top, solved half on the left (90..270 deg)
theta_plot = theta + pi/2;

ringNames = arrayfun(@(k) sprintf('r%d', k), rings, 'UniformOutput', false);

light_figure;
hold on;
axis equal;
axis off;
title(['Node temperatures [degC] (rings: ' strjoin(ringNames, ', ') ')']);

% Thin circles for orientation
tt = linspace(0, 2*pi, 400);
for ir = rings
    plot(r(ir)*cos(tt), r(ir)*sin(tt), 'k-', 'LineWidth', 0.8);
end

% Small shift so the text does not sit on top of the node marker
offset = 0.05*(rO - rI);

for g = 1:Nth
    th = theta_plot(g);

    for ir = rings
        i    = (g-1)*Nr + ir;       % index in T
        Tval = T(i);

        x = r(ir)*cos(th);
        y = r(ir)*sin(th);

        plot(x, y, 'k.', 'MarkerSize', 10);

        text(x - offset*sin(th), y + offset*cos(th), sprintf('%.1f', Tval), ...
            'Color', 'k', 'FontSize', 9, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

% Margin so the labels are not clipped
xlim(1.12*[-rO rO]);
ylim(1.12*[-rO rO]);
hold off;
end
