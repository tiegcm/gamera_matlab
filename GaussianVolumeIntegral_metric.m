function I = GaussianVolumeIntegral_metric_vec(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,...
                                               W,B1,B2,B3,B4,B5,B6,B7,B8,...
                                               D1_1,D1_2,D1_3,D1_4,D1_5,D1_6,D1_7,D1_8,...
                                               D2_1,D2_2,D2_3,D2_4,D2_5,D2_6,D2_7,D2_8,...
                                               D3_1,D3_2,D3_3,D3_4,D3_5,D3_6,D3_7,D3_8)
%%
% Gaussian-Legrande volume integral code in a general hexahedron formed by
% 8 note points (X1, X2, X3, X4, X5, X6, X7, X8)

% % Gaussian points between (-1 1)
% a =[0.1252334085  0.3678314989  0.5873179542 0.7699026741  0.9041172563  0.9815606342 -0.1252334085 -0.3678314989 -0.5873179542 -0.7699026741 -0.9041172563 -0.9815606342];
% % Gaussian weights for the corresponding quadrature points
% wt =[0.2491470458  0.2334925365  0.2031674267 0.1600783285  0.1069393259  0.0471753363 0.2491470458  0.2334925365  0.2031674267 0.1600783285  0.1069393259  0.0471753363];
% 
% N = length(a);
% 
% W  = zeros(N,N,N);
% B1 = zeros(N,N,N);
% B2 = zeros(N,N,N);
% B3 = zeros(N,N,N);
% B4 = zeros(N,N,N);
% B5 = zeros(N,N,N);
% B6 = zeros(N,N,N);
% B7 = zeros(N,N,N);
% B8 = zeros(N,N,N);
% 
% D1_1 = zeros(N,N,N);
% D1_2 = zeros(N,N,N);
% D1_3 = zeros(N,N,N);
% D1_4 = zeros(N,N,N);
% D1_5 = zeros(N,N,N);
% D1_6 = zeros(N,N,N);
% D1_7 = zeros(N,N,N);
% D1_8 = zeros(N,N,N);
% 
% D2_1 = zeros(N,N,N);
% D2_2 = zeros(N,N,N);
% D2_3 = zeros(N,N,N);
% D2_4 = zeros(N,N,N);
% D2_5 = zeros(N,N,N);
% D2_6 = zeros(N,N,N);
% D2_7 = zeros(N,N,N);
% D2_8 = zeros(N,N,N);
% 
% D3_1 = zeros(N,N,N);
% D3_2 = zeros(N,N,N);
% D3_3 = zeros(N,N,N);
% D3_4 = zeros(N,N,N);
% D3_5 = zeros(N,N,N);
% D3_6 = zeros(N,N,N);
% D3_7 = zeros(N,N,N);
% D3_8 = zeros(N,N,N);
% 
% 
% % test x,y,z coordinates - 8 points with a specific order, in the mapped
% % coordinate system {(-1 1),(-1 1),(-1 1)}, the 8 node points should follow
% % the order as: (-1,-1,-1) -> (+1,-1,-1) -> (+1,+1,-1) -> (-1,+1,-1) -> (-1,-1,+1) -> (+1,-1,+1) -> (+1,+1,+1) -> (-1,+1,+1)
% % 
% x1= -2; y1= 0; z1= 0;
% x2= 1; y2= 0; z2= 0;
% x3= 1; y3= 1; z3= 0;
% x4= 0; y4= 1; z4= 0;
% x5= 0; y5= 0; z5= 1;
% x6= 1; y6= 0; z6= 1;
% x7= 1; y7= 1; z7= 1;
% x8= 0; y8= 1; z8= 1;
% %%
%  I=0;
%  for i=1:N
%      for j=1:N
%          for k=1:N
%              
%              % weight at the (i,j,k) Gaussian Point
%              W(i,j,k) = wt(i)*wt(j)*wt(k);
%              
%              % mapped Gaussian points (P,Q,R) at i,j,k
%              B1(i,j,k)=(1-a(i))*(1-a(j))*(1-a(k)); 
%              B2(i,j,k)=(1+a(i))*(1-a(j))*(1-a(k));
%              B3(i,j,k)=(1+a(i))*(1+a(j))*(1-a(k));
%              B4(i,j,k)=(1-a(i))*(1+a(j))*(1-a(k));
%              B5(i,j,k)=(1-a(i))*(1-a(j))*(1+a(k));
%              B6(i,j,k)=(1+a(i))*(1-a(j))*(1+a(k));
%              B7(i,j,k)=(1+a(i))*(1+a(j))*(1+a(k));
%              B8(i,j,k)=(1-a(i))*(1+a(j))*(1+a(k));
%                           
%              % Jacobian determinant at the (i,j,k) Gaussian Point
%              D1_1(i,j,k)=(-1)*(1-a(j))*(1-a(k));
%              D1_2(i,j,k)=(+1)*(1-a(j))*(1-a(k));
%              D1_3(i,j,k)=(+1)*(1+a(j))*(1-a(k));
%              D1_4(i,j,k)=(-1)*(1+a(j))*(1-a(k));
%              D1_5(i,j,k)=(-1)*(1-a(j))*(1+a(k));
%              D1_6(i,j,k)=(+1)*(1-a(j))*(1+a(k));
%              D1_7(i,j,k)=(+1)*(1+a(j))*(1+a(k));
%              D1_8(i,j,k)=(-1)*(1+a(j))*(1+a(k));
%              
%              D2_1(i,j,k)=(1-a(i))*(-1)*(1-a(k));
%              D2_2(i,j,k)=(1+a(i))*(-1)*(1-a(k));
%              D2_3(i,j,k)=(1+a(i))*(+1)*(1-a(k));
%              D2_4(i,j,k)=(1-a(i))*(+1)*(1-a(k));
%              D2_5(i,j,k)=(1-a(i))*(-1)*(1+a(k));
%              D2_6(i,j,k)=(1+a(i))*(-1)*(1+a(k));
%              D2_7(i,j,k)=(1+a(i))*(+1)*(1+a(k));
%              D2_8(i,j,k)=(1-a(i))*(+1)*(1+a(k));                           
% 
%              D3_1(i,j,k)=(1-a(i))*(1-a(j))*(-1);
%              D3_2(i,j,k)=(1+a(i))*(1-a(j))*(-1);
%              D3_3(i,j,k)=(1+a(i))*(1+a(j))*(-1);
%              D3_4(i,j,k)=(1-a(i))*(1+a(j))*(-1);
%              D3_5(i,j,k)=(1-a(i))*(1-a(j))*(+1);
%              D3_6(i,j,k)=(1+a(i))*(1-a(j))*(+1);
%              D3_7(i,j,k)=(1+a(i))*(1+a(j))*(+1);
%              D3_8(i,j,k)=(1-a(i))*(1+a(j))*(+1);  
% 
%          end
%      end
%  end
% %%
%  I=0;
%  for i=1:N
%      for j=1:N
%          for k=1:N
%              
%              % weight at the (i,j,k) Gaussian Point
%              W(i,j,k) = wt(i)*wt(j)*wt(k);
%         tic    

%              % mapped Gaussian points (P,Q,R) at i,j,k
             P = 1/8*( x1*B1 + x2*B2 + x3*B3 + x4*B4 + x5*B5 + x6*B6 + x7*B7 + x8*B8 );
             Q = 1/8*( y1*B1 + y2*B2 + y3*B3 + y4*B4 + y5*B5 + y6*B6 + y7*B7 + y8*B8 );
             R = 1/8*( z1*B1 + z2*B2 + z3*B3 + z4*B4 + z5*B5 + z6*B6 + z7*B7 + z8*B8 );
             
             % Jacobian determinant at the (i,j,k) Gaussian Point
             dxd1 = 1/8*( x1*D1_1 + x2*D1_2 + x3*D1_3 + x4*D1_4 + x5*D1_5 + x6*D1_6 + x7*D1_7 + x8*D1_8 );
             dyd1 = 1/8*( y1*D1_1 + y2*D1_2 + y3*D1_3 + y4*D1_4 + y5*D1_5 + y6*D1_6 + y7*D1_7 + y8*D1_8 );
             dzd1 = 1/8*( z1*D1_1 + z2*D1_2 + z3*D1_3 + z4*D1_4 + z5*D1_5 + z6*D1_6 + z7*D1_7 + z8*D1_8 );
             
             dxd2 = 1/8*( x1*D2_1 + x2*D2_2 + x3*D2_3 + x4*D2_4 + x5*D2_5 + x6*D2_6 + x7*D2_7 + x8*D2_8 );
             dyd2 = 1/8*( y1*D2_1 + y2*D2_2 + y3*D2_3 + y4*D2_4 + y5*D2_5 + y6*D2_6 + y7*D2_7 + y8*D2_8 );
             dzd2 = 1/8*( z1*D2_1 + z2*D2_2 + z3*D2_3 + z4*D2_4 + z5*D2_5 + z6*D2_6 + z7*D2_7 + z8*D2_8 );

             dxd3 = 1/8*( x1*D3_1 + x2*D3_2 + x3*D3_3 + x4*D3_4 + x5*D3_5 + x6*D3_6 + x7*D3_7 + x8*D3_8 );
             dyd3 = 1/8*( y1*D3_1 + y2*D3_2 + y3*D3_3 + y4*D3_4 + y5*D3_5 + y6*D3_6 + y7*D3_7 + y8*D3_8 );
             dzd3 = 1/8*( z1*D3_1 + z2*D3_2 + z3*D3_3 + z4*D3_4 + z5*D3_5 + z6*D3_6 + z7*D3_7 + z8*D3_8 );
                          
%             Jacob = [dxd1 dyd1 dzd1
%                      dxd2 dyd2 dzd2
%                      dxd3 dyd3 dzd3];
            J = abs( dxd1.*dyd2.*dzd3 + dxd2.*dyd3.*dzd1 + dxd3.*dyd1.*dzd2 - ...
                     dxd3.*dyd2.*dzd1 - dxd2.*dyd1.*dzd3 - dxd1.*dyd3.*dzd2) ;
               
%             J = abs(det(Jacob));
            
%             I = I+W(i,j,k)*FUNC_TEST(P(i,j,k),Q(i,j,k),R(i,j,k))*J(i,j,k);
            I = sum(1.0.*W(:).*J(:));
%          end
%      end
%  end
%  
% toc
%          sum(I(:))
