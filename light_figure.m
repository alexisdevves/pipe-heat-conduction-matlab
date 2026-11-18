function fig = light_figure()
% LIGHT_FIGURE  Create a figure with a white background and dark text/lines,
% even when MATLAB runs with the dark theme (R2025a or later).
%
%   The plots draw in black, which would be invisible on the dark theme.

fig = figure('Color', 'w');
try
    theme(fig, 'light');        % MATLAB R2025a+ only
catch
    % older MATLAB versions have no themes: nothing to do
end
set(fig, 'Color', 'w');
end