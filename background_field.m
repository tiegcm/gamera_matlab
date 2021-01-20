function [fxi,fyi,fzi,fnormi,fsqi,fnxi,fnyi,fnzi,...
          fxj,fyj,fzj,fnormj,fsqj,fnxj,fnyj,fnzj,...
          fxk,fyk,fzk,fnormk,fsqk,fnxk,fnyk,fnzk,...
          BX0_I,BY0_I,BX0_J,BY0_J,BX0_K,BY0_K]=background_field(x,y,z,BX0,BY0,BZ0,vec_i,vec_j,vec_k)
           
% This code computes the face integrals of background magnetic fields (BX0, BY0, BZ0)
% these Gaussian quadratured face magnetic field terms (all divided by face
% area) are used in the magnetic stress calculations. The line-integrated
% magnetic fields are transformed into the edge-aligned coordinate system
% (e1,e2,n*) for electric field calculations.
      
global I J K Ip1 Jp1 Kp1 ic_act jc_act kc_act if_act jf_act kf_act
global nx_total ny_total nz_total

disp('Initializing Background Fields...');

% Initialize Gaussial surface integrals on i-face:
fxi = zeros(nx_total,ny_total-1,nz_total-1);     % int BX0*dS
fyi = zeros(nx_total,ny_total-1,nz_total-1);     % int BY0*dS
fzi = zeros(nx_total,ny_total-1,nz_total-1);     % int BZ0*dS
fnormi = zeros(nx_total,ny_total-1,nz_total-1);  % int (B0 dot n)*dS 
fsqi = zeros(nx_total,ny_total-1,nz_total-1);    % int (B0^2)*dS
fnxi = zeros(nx_total,ny_total-1,nz_total-1);    % int (B0 dot n)*BX0*dS
fnyi = zeros(nx_total,ny_total-1,nz_total-1);    % int (B0 dot n)*BY0*dS
fnzi = zeros(nx_total,ny_total-1,nz_total-1);    % int (B0 dot n)*BZ0*dS

% Initialize Gaussial surface integrals on j-face:
fxj = zeros(nx_total-1,ny_total,nz_total-1);     % int BX0*dS
fyj = zeros(nx_total-1,ny_total,nz_total-1);     % int BY0*dS
fzj = zeros(nx_total-1,ny_total,nz_total-1);     % int BZ0*dS
fnormj = zeros(nx_total-1,ny_total,nz_total-1);  % int (B0 dot n)*dS 
fsqj = zeros(nx_total-1,ny_total,nz_total-1);    % int (B0^2)*dS
fnxj = zeros(nx_total-1,ny_total,nz_total-1);    % int (B0 dot n)*BX0*dS
fnyj = zeros(nx_total-1,ny_total,nz_total-1);    % int (B0 dot n)*BY0*dS
fnzj = zeros(nx_total-1,ny_total,nz_total-1);    % int (B0 dot n)*BZ0*dS

% Initialize Gaussial surface integrals on k-face:
fxk = zeros(nx_total-1,ny_total-1,nz_total);     % int BX0*dS
fyk = zeros(nx_total-1,ny_total-1,nz_total);     % int BY0*dS
fzk = zeros(nx_total-1,ny_total-1,nz_total);     % int BZ0*dS
fnormk = zeros(nx_total-1,ny_total-1,nz_total);  % int (B0 dot n)*dS 
fsqk = zeros(nx_total-1,ny_total-1,nz_total);    % int (B0^2)*dS
fnxk = zeros(nx_total-1,ny_total-1,nz_total);    % int (B0 dot n)*BX0*dS
fnyk = zeros(nx_total-1,ny_total-1,nz_total);    % int (B0 dot n)*BY0*dS
fnzk = zeros(nx_total-1,ny_total-1,nz_total);    % int (B0 dot n)*BZ0*dS

% Initialize Gaussian line integrals on i-edge:
lbxi = zeros(nx_total-1,ny_total,nz_total);      % int BX0 dot dl
lbyi = zeros(nx_total-1,ny_total,nz_total);      % int BY0 dot dl
lbzi = zeros(nx_total-1,ny_total,nz_total);      % int BZ0 dot dl

% Initialize Gaussian line integrals on j-edge:
lbxj = zeros(nx_total,ny_total-1,nz_total);      % int BX0 dot dl
lbyj = zeros(nx_total,ny_total-1,nz_total);      % int BY0 dot dl
lbzj = zeros(nx_total,ny_total-1,nz_total);      % int BZ0 dot dl

% Initialize Gaussian line integrals on k-edge:
lbxk =zeros(nx_total,ny_total,nz_total-1);       % int BX0 dot dl
lbyk =zeros(nx_total,ny_total,nz_total-1);       % int BY0 dot dl
lbzk =zeros(nx_total,ny_total,nz_total-1);       % int BZ0 dot dl

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

% use the g2int method as in the LFM code
% Gaussian quadrature on i-faces
for i=1:Ip1(end)
    for j=1:J(end)
        for k=1:K(end)
            [fxi(i,j,k),fyi(i,j,k),fzi(i,j,k),fnormi(i,j,k),fsqi(i,j,k),fnxi(i,j,k),fnyi(i,j,k),fnzi(i,j,k)] = ...
             g2int(x(i,j,k),x(i,j+1,k),x(i,j,k+1),x(i,j+1,k+1),y(i,j,k),y(i,j+1,k),y(i,j,k+1),y(i,j+1,k+1),z(i,j,k),z(i,j+1,k),z(i,j,k+1),z(i,j+1,k+1),...
                    eta,psi,wt2,BX0,BY0,BZ0);
        end
    end
end

% Gaussian quadrature on j-faces
for i=1:I(end)
    for j=1:Jp1(end)
        for k=1:K(end)
            [fxj(i,j,k),fyj(i,j,k),fzj(i,j,k),fnormj(i,j,k),fsqj(i,j,k),fnxj(i,j,k),fnyj(i,j,k),fnzj(i,j,k)] = ...
             g2int(x(i,j,k),x(i,j,k+1),x(i+1,j,k),x(i+1,j,k+1),y(i,j,k),y(i,j,k+1),y(i+1,j,k),y(i+1,j,k+1),z(i,j,k),z(i,j,k+1),z(i+1,j,k),z(i+1,j,k+1),...
                    eta,psi,wt2,BX0,BY0,BZ0);
        end
    end
end

% Gaussian quadrature on k-faces
for i=1:I(end)
    for j=1:J(end)
        for k=1:Kp1(end)
            [fxk(i,j,k),fyk(i,j,k),fzk(i,j,k),fnormk(i,j,k),fsqk(i,j,k),fnxk(i,j,k),fnyk(i,j,k),fnzk(i,j,k)] = ...
             g2int(x(i,j,k),x(i+1,j,k),x(i,j+1,k),x(i+1,j+1,k),y(i,j,k),y(i+1,j,k),y(i,j+1,k),y(i+1,j+1,k),z(i,j,k),z(i+1,j,k),z(i,j+1,k),z(i+1,j+1,k),...
                    eta,psi,wt2,BX0,BY0,BZ0);
        end
    end
end
             
% integrate the background b field on edges for electric field calculations
% i-edge background B field
% integrate B along i edge
% the gaussian line integral gives the average function along a line, so in
% the backgound B calculation it gives an average field along a line rather
% than B*dL
for i=ic_act(1):ic_act(end)
    for j=jf_act(1):jf_act(end)
        for k=kf_act(1):kf_act(end)
            lbxi(i,j,k) =  GaussianLineIntegral(BX0,x(i+1,j,k),y(i+1,j,k),z(i+1,j,k),x(i,j,k),y(i,j,k),z(i,j,k));
            lbyi(i,j,k) =  GaussianLineIntegral(BY0,x(i+1,j,k),y(i+1,j,k),z(i+1,j,k),x(i,j,k),y(i,j,k),z(i,j,k));
            lbzi(i,j,k) =  GaussianLineIntegral(BZ0,x(i+1,j,k),y(i+1,j,k),z(i+1,j,k),x(i,j,k),y(i,j,k),z(i,j,k));
        end
    end
end
% integrate B along j edge
for i=if_act(1):if_act(end)
    for j=jc_act(1):jc_act(end)
        for k=kf_act(1):kf_act(end)
            lbxj(i,j,k) = GaussianLineIntegral(BX0,x(i,j+1,k),y(i,j+1,k),z(i,j+1,k),x(i,j,k),y(i,j,k),z(i,j,k));
            lbyj(i,j,k) = GaussianLineIntegral(BY0,x(i,j+1,k),y(i,j+1,k),z(i,j+1,k),x(i,j,k),y(i,j,k),z(i,j,k));
            lbzj(i,j,k) = GaussianLineIntegral(BZ0,x(i,j+1,k),y(i,j+1,k),z(i,j+1,k),x(i,j,k),y(i,j,k),z(i,j,k));
        end
    end
end
% integrate B along k edge
for i=if_act(1):if_act(end)
    for j=jf_act(1):jf_act(end)
        for k=kc_act(1):kc_act(end)
            lbxk(i,j,k) = GaussianLineIntegral(BX0,x(i,j,k+1),y(i,j,k+1),z(i,j,k+1),x(i,j,k),y(i,j,k),z(i,j,k));
            lbyk(i,j,k) = GaussianLineIntegral(BY0,x(i,j,k+1),y(i,j,k+1),z(i,j,k+1),x(i,j,k),y(i,j,k),z(i,j,k));
            lbzk(i,j,k) = GaussianLineIntegral(BZ0,x(i,j,k+1),y(i,j,k+1),z(i,j,k+1),x(i,j,k),y(i,j,k),z(i,j,k));
        end
    end
end

% transform the line-integrated background field vectors into the
% edge-aligned coordinate system for electric field calculation vx(b+B0)
% background b field transformed into K-edge coordinate
BX0_K(if_act,jf_act,kc_act) = lbxk(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,1) + ...
                              lbyk(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,2) + ...
                              lbzk(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,3);
BY0_K(if_act,jf_act,kc_act) = lbxk(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,4) + ...
                              lbyk(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,5) + ...
                              lbzk(if_act,jf_act,kc_act).*vec_k(if_act,jf_act,kc_act,6);
% transform background B to the I-edge system
BX0_I(ic_act,jf_act,kf_act) = lbxi(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,1) + ...
                              lbyi(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,2) + ...
                              lbzi(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,3);
BY0_I(ic_act,jf_act,kf_act) = lbxi(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,4) + ...
                              lbyi(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,5) + ...
                              lbzi(ic_act,jf_act,kf_act).*vec_i(ic_act,jf_act,kf_act,6);
% transform background B field to j-edge coordinate sytem
BX0_J(if_act,jc_act,kf_act) = lbxj(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,1) + ...
                              lbyj(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,2) + ...
                              lbzj(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,3);
BY0_J(if_act,jc_act,kf_act) = lbxj(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,4) + ...
                              lbyj(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,5) + ...
                              lbzj(if_act,jc_act,kf_act).*vec_j(if_act,jc_act,kf_act,6);                          

% END of background B fields
end

