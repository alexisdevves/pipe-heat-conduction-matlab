function T_ext = T_ext_fun(theta)
% T_EXT_FUN  Outer-surface temperature of the pipe [degC] (Dirichlet boundary condition).
%
%   theta : angle [rad], measured from the TOP (north) of the pipe, so the
%           function is symmetric about the vertical axis.
%           theta = 0 -> top (hottest), theta = pi -> bottom (coldest).
%
%   T_ext(theta) = 24.5 + 15*cos(theta) - 2.5*cos(2*theta)

T_ext = 24.5 + 15*cos(theta) - 2.5*cos(2*theta);
end
