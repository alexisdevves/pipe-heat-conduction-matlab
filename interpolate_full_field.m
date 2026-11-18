function [X, Y, Tfull, r_nodes] = interpolate_full_field(T, Nr, Nth, rI, rO, theta, iM)
% INTERPOLATE_FULL_FIELD  Smooth 360-degree temperature field for plotting.
%
%   Takes the nodal solution on the half section (0..pi) and returns a dense
%   polar field over the full circle, built by mirroring about the vertical axis.
%
%   Inputs
%     T      : solution vector (Nr*Nth x 1), node id = (g-1)*Nr + ir
%     Nr,Nth : number of radial / angular nodes
%     rI,rO  : inner / outer radius [m]
%     theta  : spoke angles on the half section, 1 x Nth, from 0 to pi
%     iM     : index of the ring on the material interface
%
%   Outputs
%     X, Y   : Cartesian coordinates of the dense grid (Nr_f x Nth_full)
%     Tfull  : interpolated temperatures on that grid
%     r_nodes: radii of the node rings (1 x Nr)
%
%   NOTE: this is a visual interpolation, not part of the numerical solution.
%   The radial interpolation is done separately in each material, so the
%   change of temperature gradient at the interface (jump in conductivity)
%   is not smoothed away.

if numel(T) ~= Nr*Nth
    error('T has %d elements but Nr*Nth = %d. Run the solver again.', numel(T), Nr*Nth);
end

dr      = (rO - rI)/(Nr - 1);
r_nodes = rI + (0:Nr-1)*dr;

% Vector -> matrix (Nr x Nth): Tmat(ir,g) = T((g-1)*Nr + ir)
Tmat = reshape(T(:), Nr, Nth);

% Dense grids
Nr_f      = 450;
Nth_full  = 900;
r_f       = linspace(rI, rO, Nr_f);
theta_full = linspace(0, 2*pi, Nth_full);

% Fold the full circle onto the solved half: theta in (pi, 2pi] -> 2pi - theta
theta_map = theta_full;
mask = theta_map > pi;
theta_map(mask) = 2*pi - theta_map(mask);

% 1) Angular interpolation on every node ring  (Nr x Nth_full)
T_ang = interp1(theta, Tmat.', theta_map, 'spline').';

% 2) Radial interpolation, separately inside material A and material B
Tfull = zeros(Nr_f, Nth_full);
inA = r_f <= r_nodes(iM);
Tfull(inA, :)  = interp1(r_nodes(1:iM),  T_ang(1:iM, :),  r_f(inA),  'pchip');
Tfull(~inA, :) = interp1(r_nodes(iM:Nr), T_ang(iM:Nr, :), r_f(~inA), 'pchip');

% Polar -> Cartesian; angles are shifted by pi/2 so that theta = 0 is at the top
theta_plot = theta_full + pi/2;
[TH, RR] = meshgrid(theta_plot, r_f);
X = RR.*cos(TH);
Y = RR.*sin(TH);
end
