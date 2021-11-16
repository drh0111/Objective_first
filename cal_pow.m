function [mpow_t, mpow_m, error, ratio] = cal_pow(Hz, Ey, mode)
%
% Description
%   Given an area (can be several layers) of Hz and Ey, this function will
%   calculate the complex power of the total field and the complex power of
%   the expected mode, and then do some analysis.
%
% Inputs
%   HZ: A two dimensional matrix.
%       It is the total Hz field, each row is a cross section (in y
%       direction) 
%
%   EY: A two dimensional matrix.
%       It is the total Ey field, each row is a cross section (in y
%       direction) 
%
%   MODE: A 2-element structure.
%       MODE.HZ: A column vector.
%           It is the eigenmode we choose to analyze the power.
%       MODE.EY: A column vector.
%           It is the eigenmode we choose to analyze the power.
%
% Outputs 
%   MPOW_A: A real value.
%       The mean value of the total complex power (real part) of each
%       layer.
%
%   MPOW_M: A real value.
%       The mean value of the chosen mode's complex power (real part) of
%       each layer.
%
%   ERROR: 2-element real vector.
%       Represents the standard derivation of the total and chosen mode's
%       complex power (of each layer) in the considered region.
%
%   RATIO: A real value.
%       It equals MPOW_M / MPOW_T.

pow_t = real(sum(conj(Hz) .* Ey, 2)); % calculate the total power for each layer
mpow_t = mean(pow_t);
error(1) = std(pow_t);

% project the total field to expected mode
m_Hz = project(Hz, mode.Hz);
m_Ey = project(Ey, mode.Ey);

pow_m = real(sum(conj(m_Hz) .* m_Ey, 2));
mpow_m = mean(pow_m);
error(2) = std(pow_m);

ratio = mpow_m / mpow_t;

fprintf(' Mean total power: %.2f\n Error: %.2f\n Mean mode power: %.2f\n Error: %.2f\n Ratio: %.2f\n ',...
    mpow_t, error(1), mpow_m, error(2), ratio);
end

function [p] = project(x, m)
% project X to the expected mode M, x is row vector, m is column vector
p = (x * conj(m)) * m.' / (m' * m);
end