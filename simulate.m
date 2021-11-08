function [Ex, Ey, Hz] = simulate(eps, thick, prob)
%
% Description
%   This function is used to simulate the optimized dielectric structure
%   using FDFD. It first expands the raw structure to include input and
%   output waveguide, then puts source to test the performance of the
%   designed dielectric structure, and draw the image to show the results.
%
% Inputs
%   EPS: 2 dimensional matrix.
%       It is the raw structure used in 'optimize' function, only contains
%       2 layers input and output waveguide.
%
%   THICK: 4-element row vector.
%       The 4 elements represent the expanded thickness in the up, down,
%       left and right direction (not pml here). 
%
%   PROB: A structure from 'form' function.
%       the information that will be extracted form PROB:
%           OMEGA: Angular frequency.
%           INPUT: Expected Hz field in the input waveguide.
%           BC: Boundary condition.
%
% Outputs
%   EX: 2-dimensional matrix, the same size as expandes EPS.
%       The simulated EX field.
%
%   EY: 2-dimensional matrix, the same size as expandes EPS.
%       The simulated Ey field.
%
%   HZ: 2-dimensional matrix, the same size as expandes EPS.
%       The simulated Hz field.


    %
    % expand EPS, get the INPUT
    %

eps = pad_eps(eps, thick);
dims = size(eps);
range_y = thick(3) + 1 : dims(2) - thick(4); % range of the input field in y
range_xo = dims(1) - ceil(thick(2)/2) : dims(1) - 1; % range of output-energy analysis
range_xi = floor(thick(1)/2) : thick(1); % range of input-energy analysis

input.Hz = zeros(1, dims(2));
input.Hz(range_y) = prob.Hz(1, :);
input.k = angle(prob.Hz(2, end/2) / prob.Hz(1, end/2));


    %
    % Simulate and analyze the results
    %
    
% simulate the model
[Ex_w, Ey_w, Hz_w] = fdfd_2d(repmat(eps(1, :), dims(1), 1), input, ...
                        prob.omega, prob.bc); % only consider the waveguide to get input power
[Ex, Ey, Hz] = fdfd_2d(eps, input, prob.omega, prob.bc); % consider the structure to get output power

% make some analysis
Hz_i = Hz_w(range_xi, range_y);
Ey_i = Ey_w(range_xi, range_y);
Hz_o = Hz(range_xo, range_y);
Ey_o = Ey(range_xo, range_y);

fprintf('Input mode analysis\n');
[in.mpow_t, in.mpow_m, in.error, in.ratio] = cal_pow(Hz_i, Ey_i, prob.in);

fprintf('\nOutput mode analysis\n');
[out.mpow_t, out.mpow_m, out.error, out.ratio] = cal_pow(Hz_o, Ey_o, prob.out);

eff = out.mpow_m / in.mpow_m;
fprintf('\nConversion effeciency: %.2f', eff);

% draw the plot
mul_plot(dims, {'\epsilon', eps}, {'|Hz|', abs(Hz)}, {'Re(Hz)', real(Hz)});

end