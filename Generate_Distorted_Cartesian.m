function [x,y,z]=Generate_Distorted_Cartesian(nx,ny,nz,v0,NO2)

% This code generates a distorted Cartesian grid betwee [0 1]x[0 1]
% INPUT: nx,ny,nz - # of active cells
%        NO2      - # of ghost cells
%        v0       - level of distortion, if v0 = 0 then no distortion

% grid spacing based on the # of active cells
deltax = 2/nx;
deltay = 2/ny;
deltaz = 2/nz;

% 1-D grids in each dimension
x1 = -1-NO2*deltax:deltax:1+NO2*deltax;
y1 = -1-NO2*deltay:deltay:1+NO2*deltay;
z1 = -1-NO2*deltaz:deltaz:1+NO2*deltaz;

% x,y,z are the cell corner grids (including ghost cells)
x=zeros(nx+NO2*2+1,ny+NO2*2+1,nz+NO2*2+1);
y=x;
z=x;

% fill in the 3-D grid
for i=1:length(x1)
    for j=1:length(y1)
        for k=1:length(z1)
            x(i,j,k) = 1.*(x1(i)+v0.*sin(pi*x1(i)).*sin(pi.*y1(j)));
            y(i,j,k) = 1.*(y1(j)+v0.*sin(pi*x1(i)).*sin(pi.*y1(j)));
            z(i,j,k) = z1(k);
        end
    end
end

% scale the grid from [-1 1]x[-1 1] to [0 1]x[0 1]
x = (x+1)/2;
y = (y+1)/2;
z = (z+1)/2;