function leftTemp = BC_external_exp(~,state)
input_parameter_file='inputs/1D_input_matrix.mat';

load(input_parameter_file,'diffusivity','time_shift','peak_temp')
%boundaryFileHeatedBlock Temperature boundary conditions for heated block example
% Temperature boundary condition is defined on the left edge of the block
% in the heated block example.
%
% loc   - application region struct passed in for information purposes
% state - solution state struct passed in for information purposes
% Copyright 2014-2016 The MathWorks, Inc.
% The temperature returned depends on the solution time.
if(isnan(state.time))
  % Returning a NaN for any component of q, g, h, r when time=NaN
  % tells the solver that the boundary conditions are functions of time.
  % The PDE Toolbox documentation discusses this requirement in more detail.
  leftTemp = NaN;
elseif(state.time <= 1000)
  %****This is where I want to modify the existing function****
  
  %'State' is originally the only input it takes and then it extracts time.
  %I want to be able to use p and polyval to easily import my fitted
  %polynomials into here instead of copy and pasting what you see below.
  
  t=state.time;
  %original values commented out, new values rewritten as avars below
  
  peak=peak_temp;%35000; %(~35000)
  a=diffusivity; %3.6e7; % "diffusivity" (4.0e7)
  b=time_shift; %1.7e-8; %time shift (1.8e-8)
  
  leftTemp=peak*(1+2*(-exp(-pi^2*a*(t-b)/0.61^2)+exp(-4*pi^2*a*(t-b)/0.61^2)-exp(-9*pi^2*a*(t-b)/0.61^2)));
else
  leftTemp = 0;
end
end