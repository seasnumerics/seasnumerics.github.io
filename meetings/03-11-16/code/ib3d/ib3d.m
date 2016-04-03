% ib3D.m
% This script is the main program.

global dt nv N h rho mu ip im a l1 l2 K;
initialize_3d
init_a
timestep = 400;
if ~exist('op_sphere', 'dir')
  mkdir('op_sphere');
end
folder_path = 'op_sphere/';

for k = 1 : timestep
  disp(['Simulation timestep = ', sprintf('%03d',k-1)]);
  
  % Save Structure
  file_name = strcat(folder_path, sprintf('%03d',k-1), '.vtk');
  vtk_cell(nv, ntmax, X, v, file_name);
  
  % Move the structure
  XX = X + .5*dt*interp(u,X);
  
  % Spead Lagrangian force to Eulerian coordinates
  ff = spread(Force(XX),XX);
  
  % Solve Navier-Stokes
  [u,uu] = fluid(u,ff);
  
  % Move the structure based on the interpolated velocity
  X = X + dt*interp(uu,XX);
end