function plot_heatmap_click(T, Nr, Nth, rI, rO, theta, iM)
% PLOT_HEATMAP_CLICK  Full 360-degree polar heatmap; click on a node to see
% its id, its (ring, spoke) indices and its temperature (data cursor).
%
%   T      : solution vector (Nr*Nth x 1)
%   Nr,Nth : number of radial / angular nodes
%   rI,rO  : inner / outer radius [m]
%   theta  : spoke angles on the half section (1 x Nth, 0..pi)
%   iM     : index of the ring on the material interface
%
%   Nodes of the mirrored (right) half carry the same id and temperature as
%   their counterparts on the solved (left) half.

[X, Y, Tfull, r_nodes] = interpolate_full_field(T, Nr, Nth, rI, rO, theta, iM);
T = T(:);

theta_half_plot = theta + pi/2;     % 90..270 deg: solved half, on the left

%% Heatmap
fig = light_figure;
hp = pcolor(X, Y, Tfull);
shading interp;
axis equal;
axis off;
colorbar;
title('Click on a node to display its temperature');
hold on;

xlim(1.05*[-rO rO]);
ylim(1.05*[-rO rO]);

set(hp, 'HitTest', 'off');          % only the nodes are clickable

%% Node lists (solved half + mirrored half), all as column vectors
ids   = reshape(1:Nr*Nth, Nr, Nth);                 % ids(ir,g)
IR    = repmat((1:Nr)', 1, Nth);                    % ring index of each node
G     = repmat(1:Nth,   Nr, 1);                     % spoke index of each node
R     = repmat(r_nodes(:), 1, Nth);                 % radius of each node
THp   = repmat(theta_half_plot, Nr, 1);             % plot angle of each node

% Solved (left) half: all nodes
x_left = R.*cos(THp);
y_left = R.*sin(THp);

% Mirrored (right) half: spokes 2..Nth-1 (the axis spokes already exist)
sel     = G >= 2 & G <= Nth-1;
x_right = R(sel).*cos(pi - THp(sel));
y_right = R(sel).*sin(pi - THp(sel));

x_nodes  = [x_left(:); x_right(:)];
y_nodes  = [y_left(:); y_right(:)];
id_nodes = [ids(:);    ids(sel)];
ir_nodes = [IR(:);     IR(sel)];
g_nodes  = [G(:);      G(sel)];
T_nodes  = T(id_nodes);

%% Draw nodes
hNodes = plot(x_nodes, y_nodes, 'k.', 'MarkerSize', 12);
uistack(hNodes, 'top');

hNodes.UserData.x  = x_nodes;
hNodes.UserData.y  = y_nodes;
hNodes.UserData.T  = T_nodes;
hNodes.UserData.id = id_nodes;
hNodes.UserData.ir = ir_nodes;
hNodes.UserData.g  = g_nodes;

%% Data cursor
dcm = datacursormode(fig);
set(dcm, 'Enable', 'on');
set(dcm, 'UpdateFcn', @(~, evt) nodeClickText(evt, hNodes));

hold off;
end

%% ===============================================================
% Callback: text shown for the node closest to the click
% ===============================================================
function txt = nodeClickText(evt, hNodes)

pos = evt.Position;
xq = pos(1);
yq = pos(2);

% Nearest node
[~, idx] = min((hNodes.UserData.x - xq).^2 + (hNodes.UserData.y - yq).^2);

txt = {
    sprintf('Node ID: %d', hNodes.UserData.id(idx))
    sprintf('(ir, g) = (%d, %d)', hNodes.UserData.ir(idx), hNodes.UserData.g(idx))
    sprintf('T = %.4f degC', hNodes.UserData.T(idx))
    };
end
