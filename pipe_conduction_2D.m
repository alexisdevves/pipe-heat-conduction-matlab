%% ===============================================================
%  2D steady-state heat conduction in a composite cylindrical wall
%  Finite differences on a polar grid, solved with Gauss-Seidel.
%
%  Only the half section 0 <= theta <= pi is solved (symmetry).
%  theta is measured from the TOP (north) of the pipe, so the
%  symmetry axis is the vertical axis.
%
%  Node numbering: i = (g-1)*Nr + ir, where
%     ir = 1..Nr   radial index (1 = inner surface, Nr = outer surface)
%     g  = 1..Nth  angular index (spoke)
%  Radial neighbours of node i are i-1 / i+1,
%  angular neighbours are i-Nr / i+Nr.
% ===============================================================
clear; close all; clc;

%% ---------------- PROBLEM DATA ----------------
% Geometry
rI = 25e-3;      % [m] inner radius
rM = 28e-3;      % [m] material interface radius (must coincide with a grid ring)
rO = 40e-3;      % [m] outer radius

% Grid
Nr  = 6;         % number of radial nodes (rings), Nr >= 3
Nth = 19;        % number of angular nodes (spokes) on the half section

% Thermal properties
kA = 15.0;       % [W/(m K)] conductivity of material A (inner, rI..rM)
kB = 0.2;        % [W/(m K)] conductivity of material B (outer, rM..rO)
hf = 4000.0;     % [W/(m^2 K)] inner convection coefficient

% Temperatures
Tf = 97.0;       % [degC] fluid temperature (inside the pipe)
% Outer boundary: Dirichlet, T_ext(theta) -> see T_ext_fun.m

% Iteration parameters
tolerance      = 1e-9;
max_iterations = 200000;

%% ---------------- GRID ----------------
N     = Nr*Nth;                 % total number of nodes
dr    = (rO - rI)/(Nr - 1);     % radial step
dth   = pi/(Nth - 1);           % angular step
theta = linspace(0, pi, Nth);   % angles of the spokes (1 x Nth)
r     = rI + (0:Nr-1)*dr;       % radii of the rings   (1 x Nr)

% Ring that lies on the material interface
iM = round((rM - rI)/dr) + 1;
if iM < 2 || iM > Nr - 1 || abs(r(iM) - rM) > 1e-6*dr
    error(['The interface radius rM = %.4g m does not coincide with a grid ring ' ...
           '(or lies on the inner/outer surface). Adjust rM or Nr.'], rM);
end

% Draw the grid and the node numbering
plot_mesh(rI, rO, Nr, Nth);

%% ---------------- INITIAL CONDITIONS ----------------
T_ext_top    = T_ext_fun(0);       % theta = 0  (top)
T_ext_bottom = T_ext_fun(pi);      % theta = pi (bottom)

% Initial guess: mean of the fluid, top and bottom outer temperatures
T_init = (Tf + T_ext_top + T_ext_bottom)/3;
T = T_init*ones(N, 1);

% Dirichlet condition on the outer ring (never updated by the iterations)
outerNodes = Nr:Nr:N;
T(outerNodes) = T_ext_fun(theta);

fprintf('T_ext(0)=%.3f degC, T_ext(pi)=%.3f degC\n', T_ext_top, T_ext_bottom);
fprintf('T_init=%.3f degC\n', T_init);

%% ---------------- DISCRETIZATION COEFFICIENTS ----------------
% Per ring (1 x Nr)
aN = (r + dr/2)./(r*dr^2);      % radial, outwards
aS = (r - dr/2)./(r*dr^2);      % radial, inwards
aW = 1./(r.^2*dth^2);           % angular (theta-)
aE = aW;                        % angular (theta+)
aP = aN + aS + aW + aE;         % central coefficient of a generic node

a_conv = hf/(kA*dr);            % convection term at the inner ring

%% ---------------- GAUSS-SEIDEL (in place) ----------------
fprintf('\n--- Gauss-Seidel iteration ---\n');

converged = false;
for iter = 1:max_iterations
    T_old = T;                          % for the convergence check

    for ir = 1:Nr-1                     % rings (outer ring Nr is Dirichlet)
        for g = 1:Nth                   % spokes

            i = (g-1)*Nr + ir;          % global node index

            % Angular neighbours; at theta = 0 and theta = pi the missing
            % neighbour is replaced by its mirror image (ghost node,
            % dT/dtheta = 0).
            if g == 1
                T_m = T(i + Nr);        % ghost for i-Nr
                T_p = T(i + Nr);
            elseif g == Nth
                T_m = T(i - Nr);
                T_p = T(i - Nr);        % ghost for i+Nr
            else
                T_m = T(i - Nr);
                T_p = T(i + Nr);
            end

            if ir == 1
                % ---- inner ring: convection with the fluid ----
                T_out = T(i + 1);
                AP  = a_conv + aW(ir) + aN(ir);
                RHS = a_conv*Tf + 0.5*aW(ir)*(T_p + T_m) + aN(ir)*T_out;

            elseif ir == iM
                % ---- material interface (A inside, B outside) ----
                T_in  = T(i - 1);
                T_out = T(i + 1);
                AP  = kB*aN(ir) + kA*aS(ir) + (kA + kB)*aW(ir);
                RHS = kB*aN(ir)*T_out + kA*aS(ir)*T_in ...
                    + 0.5*(kA + kB)*aW(ir)*(T_m + T_p);

            else
                % ---- generic interior node (conductivity cancels) ----
                T_in  = T(i - 1);
                T_out = T(i + 1);
                AP  = aP(ir);
                RHS = aW(ir)*T_p + aE(ir)*T_m + aN(ir)*T_out + aS(ir)*T_in;
            end

            T(i) = RHS/AP;
        end
    end

    % ---- convergence check ----
    err = max(abs(T - T_old));
    if err < tolerance
        converged = true;
        break;
    end
end

if ~converged
    error(['Gauss-Seidel did NOT converge in %d iterations (err = %.3e). ' ...
           'Results are not final, so no plots are produced.'], ...
           max_iterations, err);
end
fprintf('Gauss-Seidel converged in %d iterations (err = %.3e)\n', iter, err);

%% ---------------- RESULTS ----------------
fprintf('\nFinal temperature vector T (%dx1):\n', N);
disp(T);
fprintf('Tmin = %.4f degC, Tmax = %.4f degC\n', min(T), max(T));

% Rings shown in the node-temperature plot: inner, interface, middle, outer
rings = unique([1, iM, round((iM + Nr)/2), Nr]);

plot_heatmap_full(T, Nr, Nth, rI, rO, theta, iM);
plot_node_temperatures(T, Nr, Nth, rI, rO, theta, rings);
plot_heatmap_click(T, Nr, Nth, rI, rO, theta, iM);
