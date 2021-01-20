function [xc,yc,zc,...
          xi,yi,zi,...
          xj,yj,zj,...
          xk,yk,zk,...
          dx,dy,dz,...
          di,dj,dk,...
          i_face_area,j_face_area,k_face_area,volume]=metric_grid(x,y,z)

% This code computes the grid locations (center, face), face area and cell
% volume, the latter two using 12-th order Gaussian quadrature methods for
% more accurate metric calculations
%
% INPUT: x,y,z - cell cornors
  
global I J K Ip1 Jp1 Kp1

disp('Computing Cell-Center Locations...');

% locations for cell centers
xc(I,J,K) = 0.125*( x(I,J,K) + x(I,J,K+1) + x(I+1,J,K+1) + x(I+1,J,K) + x(I,J+1,K) + x(I,J+1,K+1) + x(I+1,J+1,K+1) + x(I+1,J+1,K) );
yc(I,J,K) = 0.125*( y(I,J,K) + y(I,J,K+1) + y(I+1,J,K+1) + y(I+1,J,K) + y(I,J+1,K) + y(I,J+1,K+1) + y(I+1,J+1,K+1) + y(I+1,J+1,K) );
zc(I,J,K) = 0.125*( z(I,J,K) + z(I,J,K+1) + z(I+1,J,K+1) + z(I+1,J,K) + z(I,J+1,K) + z(I,J+1,K+1) + z(I+1,J+1,K+1) + z(I+1,J+1,K) );
% dx, dy, dz in each cell
dx(I,J,K) = x(I+1,J,K) - x(I,J,K);
dy(I,J,K) = y(I,J+1,K) - y(I,J,K);
dz(I,J,K) = z(I,J,K+1) - z(I,J,K);

%~~~~~~~~~~~~~~~~~~~~Initial Face Center Calculations~~~~~~~~~~~~~~~~~~~~~~
% John's method: simply using four-point averages of corner grids

% locations for i-face centers -  where bi is defined - this is a simple
% calculation using four-point averaged locations of cell corners. The
% values of xi,yi,zi will be replaced by the Gaussian quadrature next. This
% is the old calculation in LFM, here it severs as initialization for x{ijk}
yi(Ip1,J,K) = 0.25*( y(Ip1,J,K) + y(Ip1,J+1,K) + y(Ip1,J,K+1) + y(Ip1,J+1,K+1) );
xi(Ip1,J,K) = 0.25*( x(Ip1,J,K) + x(Ip1,J+1,K) + x(Ip1,J,K+1) + x(Ip1,J+1,K+1) );
zi(Ip1,J,K) = 0.25*( z(Ip1,J,K) + z(Ip1,J+1,K) + z(Ip1,J,K+1) + z(Ip1,J+1,K+1) );

% locations for j-face centers - where bj is defined 
% calculation using four-point averaged locations of cell corners. The
% values of xj,yj,zj will be replaced by the Gaussian quadrature next.
xj(I,Jp1,K) = 0.25*( x(I,Jp1,K) + x(I+1,Jp1,K) + x(I,Jp1,K+1) + x(I+1,Jp1,K+1));
yj(I,Jp1,K) = 0.25*( y(I,Jp1,K) + y(I+1,Jp1,K) + y(I,Jp1,K+1) + y(I+1,Jp1,K+1));
zj(I,Jp1,K) = 0.25*( z(I,Jp1,K) + z(I+1,Jp1,K) + z(I,Jp1,K+1) + z(I+1,Jp1,K+1));

% locations for k-face centers - where bk is defined 
% calculation using four-point averaged locations of cell corners. The
% values of xk,yk,zk will be replaced by the Gaussian quadrature next.
xk(I,J,Kp1) = 0.25*( x(I,J,Kp1) + x(I+1,J,Kp1) + x(I,J+1,Kp1) + x(I+1,J+1,Kp1));
yk(I,J,Kp1) = 0.25*( y(I,J,Kp1) + y(I+1,J,Kp1) + y(I,J+1,Kp1) + y(I+1,J+1,Kp1));
zk(I,J,Kp1) = 0.25*( z(I,J,Kp1) + z(I+1,J,Kp1) + z(I,J+1,Kp1) + z(I+1,J+1,Kp1));

%~~~~~~~~~~~~~~~~~~~Accurate Face Center Calculations~~~~~~~~~~~~~~~~~~~~~~
% Kareem's method: using high-order gaussian integrals to calculate the
% face center locations more accurately especially when cells are very
% distorted near the spherical pole axis

% Gaussian positions for 2-D Gaussian Quadrature (12x12 point)
a =[0.1252334085  0.3678314989  0.5873179542
    0.7699026741  0.9041172563  0.9815606342
    -0.1252334085 -0.3678314989 -0.5873179542
    -0.7699026741 -0.9041172563 -0.9815606342];
% Gaussian weights for 2-D Gaussian Quadrature (12x12 point)
wt =[0.2491470458  0.2334925365  0.2031674267
    0.1600783285  0.1069393259  0.0471753363
    0.2491470458  0.2334925365  0.2031674267
    0.1600783285  0.1069393259  0.0471753363];
% variables for 2-D Gaussian quadrature
eta = zeros(12,12);
psi = eta;
wt2 = eta;
for i=1:12
    for j=1:12
        eta(i,j) = 0.5.*(1. + a(i));
        psi(i,j) = 0.5.*(1. + a(j));
        wt2(i,j) = wt(i).*wt(j);
    end
end

% initialize variables for face-area (face centered) and volume (cell centered)
i_face_area = xi*0;
j_face_area = xj*0;
k_face_area = xk*0;
volume = xc*0;

disp('Computing Face-Center Locations and Area...');
% compute i-face areas (i_face_area) and i-face centers (xi,yi,zi)
% note the index difference in the three dimensions
for i=1:Ip1(end)
    for j=1:J(end)
        for k=1:K(end)
            [i_face_area(i,j,k),xi(i,j,k),yi(i,j,k),zi(i,j,k)] = ...
                GaussianFaceIntegral_metric(x(i,j,k),x(i,j+1,k),x(i,j,k+1),x(i,j+1,k+1),y(i,j,k),y(i,j+1,k),y(i,j,k+1),y(i,j+1,k+1),z(i,j,k),z(i,j+1,k),z(i,j,k+1),z(i,j+1,k+1),eta,psi,wt2);
        end
    end
end
% compute j-face areas (j_face_area) and j-face centers (xj,yj,zj)
% note the index difference in the three dimensions
for i=1:I(end)
    for j=1:Jp1(end)
        for k=1:K(end)
            [j_face_area(i,j,k),xj(i,j,k),yj(i,j,k),zj(i,j,k)] = ...
                GaussianFaceIntegral_metric(x(i,j,k),x(i,j,k+1),x(i+1,j,k),x(i+1,j,k+1),y(i,j,k),y(i,j,k+1),y(i+1,j,k),y(i+1,j,k+1),z(i,j,k),z(i,j,k+1),z(i+1,j,k),z(i+1,j,k+1),eta,psi,wt2);
        end
    end
end
% compute k-face areas (k_face_area) and k-face centers (xk,yk,zk)
% note the index difference in the three dimensions
for i=1:I(end)
    for j=1:J(end)
        for k=1:Kp1(end)
            [k_face_area(i,j,k),xk(i,j,k),yk(i,j,k),zk(i,j,k)] = ...
                GaussianFaceIntegral_metric(x(i,j,k),x(i+1,j,k),x(i,j+1,k),x(i+1,j+1,k),y(i,j,k),y(i+1,j,k),y(i,j+1,k),y(i+1,j+1,k),z(i,j,k),z(i+1,j,k),z(i,j+1,k),z(i+1,j+1,k),eta,psi,wt2);
        end
    end
end

disp('Computing Volume...')
% calculate the volume use a 12-th order 3-D Gaussian quadrature
% Gaussian points between (-1 1)
a =[0.1252334085  0.3678314989  0.5873179542 0.7699026741  0.9041172563  0.9815606342 -0.1252334085 -0.3678314989 -0.5873179542 -0.7699026741 -0.9041172563 -0.9815606342];
% Gaussian weights for the corresponding quadrature points
wt =[0.2491470458  0.2334925365  0.2031674267 0.1600783285  0.1069393259  0.0471753363 0.2491470458  0.2334925365  0.2031674267 0.1600783285  0.1069393259  0.0471753363];

N = length(a);

W  = zeros(N,N,N);
B1 = zeros(N,N,N);
B2 = zeros(N,N,N);
B3 = zeros(N,N,N);
B4 = zeros(N,N,N);
B5 = zeros(N,N,N);
B6 = zeros(N,N,N);
B7 = zeros(N,N,N);
B8 = zeros(N,N,N);

D1_1 = zeros(N,N,N);
D1_2 = zeros(N,N,N);
D1_3 = zeros(N,N,N);
D1_4 = zeros(N,N,N);
D1_5 = zeros(N,N,N);
D1_6 = zeros(N,N,N);
D1_7 = zeros(N,N,N);
D1_8 = zeros(N,N,N);

D2_1 = zeros(N,N,N);
D2_2 = zeros(N,N,N);
D2_3 = zeros(N,N,N);
D2_4 = zeros(N,N,N);
D2_5 = zeros(N,N,N);
D2_6 = zeros(N,N,N);
D2_7 = zeros(N,N,N);
D2_8 = zeros(N,N,N);

D3_1 = zeros(N,N,N);
D3_2 = zeros(N,N,N);
D3_3 = zeros(N,N,N);
D3_4 = zeros(N,N,N);
D3_5 = zeros(N,N,N);
D3_6 = zeros(N,N,N);
D3_7 = zeros(N,N,N);
D3_8 = zeros(N,N,N);

% integrate over all finite-volume (controlled) cells
for i=1:N
    for j=1:N
        for k=1:N
            
            % weight at the (i,j,k) Gaussian Point
            W(i,j,k) = wt(i)*wt(j)*wt(k);
            
            % mapped Gaussian points (P,Q,R) at i,j,k
            B1(i,j,k)=(1-a(i))*(1-a(j))*(1-a(k));
            B2(i,j,k)=(1+a(i))*(1-a(j))*(1-a(k));
            B3(i,j,k)=(1+a(i))*(1+a(j))*(1-a(k));
            B4(i,j,k)=(1-a(i))*(1+a(j))*(1-a(k));
            B5(i,j,k)=(1-a(i))*(1-a(j))*(1+a(k));
            B6(i,j,k)=(1+a(i))*(1-a(j))*(1+a(k));
            B7(i,j,k)=(1+a(i))*(1+a(j))*(1+a(k));
            B8(i,j,k)=(1-a(i))*(1+a(j))*(1+a(k));
            
            % Jacobian determinant at the (i,j,k) Gaussian Point
            D1_1(i,j,k)=(-1)*(1-a(j))*(1-a(k));
            D1_2(i,j,k)=(+1)*(1-a(j))*(1-a(k));
            D1_3(i,j,k)=(+1)*(1+a(j))*(1-a(k));
            D1_4(i,j,k)=(-1)*(1+a(j))*(1-a(k));
            D1_5(i,j,k)=(-1)*(1-a(j))*(1+a(k));
            D1_6(i,j,k)=(+1)*(1-a(j))*(1+a(k));
            D1_7(i,j,k)=(+1)*(1+a(j))*(1+a(k));
            D1_8(i,j,k)=(-1)*(1+a(j))*(1+a(k));
            
            D2_1(i,j,k)=(1-a(i))*(-1)*(1-a(k));
            D2_2(i,j,k)=(1+a(i))*(-1)*(1-a(k));
            D2_3(i,j,k)=(1+a(i))*(+1)*(1-a(k));
            D2_4(i,j,k)=(1-a(i))*(+1)*(1-a(k));
            D2_5(i,j,k)=(1-a(i))*(-1)*(1+a(k));
            D2_6(i,j,k)=(1+a(i))*(-1)*(1+a(k));
            D2_7(i,j,k)=(1+a(i))*(+1)*(1+a(k));
            D2_8(i,j,k)=(1-a(i))*(+1)*(1+a(k));
            
            D3_1(i,j,k)=(1-a(i))*(1-a(j))*(-1);
            D3_2(i,j,k)=(1+a(i))*(1-a(j))*(-1);
            D3_3(i,j,k)=(1+a(i))*(1+a(j))*(-1);
            D3_4(i,j,k)=(1-a(i))*(1+a(j))*(-1);
            D3_5(i,j,k)=(1-a(i))*(1-a(j))*(+1);
            D3_6(i,j,k)=(1+a(i))*(1-a(j))*(+1);
            D3_7(i,j,k)=(1+a(i))*(1+a(j))*(+1);
            D3_8(i,j,k)=(1-a(i))*(1+a(j))*(+1);
            
        end
    end
end

for i=1:I(end)
    for j=1:J(end)
        for k=1:K(end)
            
            volume(i,j,k) = GaussianVolumeIntegral_metric(x(i,j,k),y(i,j,k),z(i,j,k),x(i+1,j,k),y(i+1,j,k),z(i+1,j,k)...
                                                         ,x(i+1,j+1,k),y(i+1,j+1,k),z(i+1,j+1,k),x(i,j+1,k),y(i,j+1,k),z(i,j+1,k)...
                                                         ,x(i,j,k+1),y(i,j,k+1),z(i,j,k+1),x(i+1,j,k+1),y(i+1,j,k+1),z(i+1,j,k+1)...
                                                         ,x(i+1,j+1,k+1),y(i+1,j+1,k+1),z(i+1,j+1,k+1),x(i,j+1,k+1),y(i,j+1,k+1),z(i,j+1,k+1),...
                                                          W,B1,B2,B3,B4,B5,B6,B7,B8,...
                                                          D1_1,D1_2,D1_3,D1_4,D1_5,D1_6,D1_7,D1_8,...
                                                          D2_1,D2_2,D2_3,D2_4,D2_5,D2_6,D2_7,D2_8,...
                                                          D3_1,D3_2,D3_3,D3_4,D3_5,D3_6,D3_7,D3_8);
            
        end
    end
end

% average i-unit vector in cell (I,J,K)
di_x(I,J,K) = xi(I+1,J,K) - xi(I,J,K);
di_y(I,J,K) = yi(I+1,J,K) - yi(I,J,K);
di_z(I,J,K) = zi(I+1,J,K) - zi(I,J,K);

% average j-unit vector in cell (I,J,K)
dj_x(I,J,K) = xj(I,J+1,K) - xj(I,J,K);
dj_y(I,J,K) = yj(I,J+1,K) - yj(I,J,K);
dj_z(I,J,K) = zj(I,J+1,K) - zj(I,J,K);

% average k-unit vector in cell (I,J,K)
dk_x(I,J,K) = xk(I,J,K+1) - xk(I,J,K);
dk_y(I,J,K) = yk(I,J,K+1) - yk(I,J,K);
dk_z(I,J,K) = zk(I,J,K+1) - zk(I,J,K);

% estimation of mean cell length in the i, j, k directions for CFL calculations
di(I,J,K)=sqrt(di_x(I,J,K).^2+di_y(I,J,K).^2+di_z(I,J,K).^2);
dj(I,J,K)=sqrt(dj_x(I,J,K).^2+dj_y(I,J,K).^2+dj_z(I,J,K).^2);
dk(I,J,K)=sqrt(dk_x(I,J,K).^2+dk_y(I,J,K).^2+dk_z(I,J,K).^2);

end