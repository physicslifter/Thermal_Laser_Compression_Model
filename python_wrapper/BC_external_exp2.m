function leftTemp = BC_external_exp2(~,state)
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
  
  %This is where we define our boundary conditions or "initial input". 
  %right now, its a square wave so peak temp is only free parameter
  peak=36000; 
  start_time=2e-08; 
  
  t=state.time; %print time to see live progress
  if t < 4e-8 && t > start_time
      leftTemp = peak;
  else
      leftTemp=2000;
  end    

else
  leftTemp = 2000;
end
end