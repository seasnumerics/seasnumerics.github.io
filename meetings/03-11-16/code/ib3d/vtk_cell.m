function vtk_cell(nv, nt, x, v, filename)
% vtk_cell saves a 3-D array into VTK format.
%  vtk_cell(nv, nt, x, v, filename) saves a 3-D array(xx) of the
%  size (N,3) into 'filename.vtk' format.
%
% Syntax: vtk_cell(nv, nt, x, v, filename)
%
% Inputs:
%    nv: number of vertices.
%    nt: number of triangles.
%    x: vertex coordinates.
%    v: which three vertex to form a triangle.
%    filename: string filename.
%
% Outputs:
%    None
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Author: Chen-Hung Wu (chw@mit.edu)
% Mar 2016; Last revision: 19-Mar-2016

%------------- BEGIN CODE --------------

% Setting 0 for the amount small enough
x(abs(x)<1e-16) = 0;

% Open and write a VTK file
fid = fopen(filename, 'wt');
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK output\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, '\n');

% Dataset : unstructured_Grid
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

% Input the data points
fprintf(fid, 'POINTS   %d  float', nv);
fprintf(fid, '\n');
for k=1:nv
    fprintf(fid, '%d  ', x(k,:));
    fprintf(fid, '\n');
end

% Input the CELLS : Triangles among three points
fprintf(fid, '\n');
fprintf(fid, 'CELLS   %d\t%d\n', nt,4*nt);
for k=1:nt
    fprintf(fid, '3\t');
    fprintf(fid, '%d  ', v(k,:)-1); % -1 because the data in VTK starts from 0
    fprintf(fid, '\n');
end

% Input the CELLS_TYPES
fprintf(fid, '\n');
fprintf(fid, 'CELL_TYPES   %d\n', nt);
for k=1:nt
    fprintf(fid, '5\n');
end

fclose(fid);

%------------- END OF CODE --------------
end