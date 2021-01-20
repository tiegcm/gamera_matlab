function [face_area,face_center_x,face_center_y,face_center_z] = GaussianFaceIntegral_metric(x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3,eta,psi,wt2)
%%

dx(1) = x1-x0;
dx(2) = x2-x0;
dx(3) = x3+x0-x2-x1;
dy(1) = y1-y0;
dy(2) = y2-y0;
dy(3) = y3+y0-y2-y1;
dz(1) = z1-z0;
dz(2) = z2-z0;
dz(3) = z3+z0-z2-z1;

x = x0 + dx(1).*eta + dx(2).*psi + dx(3).*eta.*psi;
y = y0 + dy(1).*eta + dy(2).*psi + dy(3).*eta.*psi;
z = z0 + dz(1).*eta + dz(2).*psi + dz(3).*eta.*psi;

xeta = dx(1) + dx(3).*psi;
yeta = dy(1) + dy(3).*psi;
zeta = dz(1) + dz(3).*psi;

xpsi = dx(2) + dx(3).*eta;
ypsi = dy(2) + dy(3).*eta;
zpsi = dz(2) + dz(3).*eta;

xn = 1.*(yeta.*zpsi-zeta.*ypsi);
yn = 1.*(zeta.*xpsi-xeta.*zpsi);
zn = 1.*(xeta.*ypsi-xpsi.*yeta);

% bx = BX0(x,y,z);
% by = BY0(x,y,z);
% bz = BZ0(x,y,z);

rnsq = xn.*xn+yn.*yn+zn.*zn;
rn = sqrt(rnsq);
area = .25.*rn.*wt2;

face_area = sum(area(:).*1.0);
face_center_x = sum(area(:).*x(:))./face_area;
face_center_y = sum(area(:).*y(:))./face_area;
face_center_z = sum(area(:).*z(:))./face_area;

% fx = sum(area(:).*bx(:));
% fy = sum(area(:).*by(:));
% fz = sum(area(:).*bz(:));
% 
% fx2 = sum(area(:).*bx(:).*bx(:));
% fy2 = sum(area(:).*by(:).*by(:));
% fz2 = sum(area(:).*bz(:).*bz(:));
% 
% fxy = sum(area(:).*bx(:).*by(:));
% fyz = sum(area(:).*by(:).*bz(:));
% fzx = sum(area(:).*bz(:).*bx(:));

end