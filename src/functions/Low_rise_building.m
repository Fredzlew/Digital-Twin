% HEYYYYY
%clc; clear; 
close all;
%% Parameters

% Dimensions of the building
H = 3.5; % Story heigth (m)
W = 8; % Width of building (m)
L = 9; % Depth of building (m) 
w = 0.5; % Width of column sections at each storey level (m)
d = 0.5; % Depth of column sections at each storey level (m)
h = 0.8; % Height of beam sections at each storey level (m)
t = 0.35; % Width of beam sections at each storey level (m)
rho_c = 2400; % Mass density of concrete (kg/m3)
rho_g = 2600; % Mass density of glas (kg/m3)

% %Jakob og Rachels tal
% H = 2.5; % Story heigth (m)
% W = 5; % Width of building (m)
% L = 7.5; % Depth of building (m) 
% w = 0.3; % Width of column sections at each storey level (m)
% d = 0.6; % Depth of column sections at each storey level (m)
% h = 0.5; % Height of beam sections at each storey level (m)
% t = 0.25; % Width of beam sections at each storey level (m)
% rho_c = 2400; % Mass density of concrete (kg/m3)
% rho_g = 2600; % Mass density of glas (kg/m3)

% Hjælpelærers tal
% H = 3.5; % Story heigth (m)
% W = 4; % Width of building (m)
% L = 6; % Depth of building (m) 
% w = 0.3; % Width of column sections at each storey level (m)
% d = 0.5; % Depth of column sections at each storey level (m)
% h = 0.5; % Height of beam sections at each storey level (m)
% t = 0.25; % Width of beam sections at each storey level (m)
% rho_c = 2400; % Mass density of concrete (kg/m3)
% rho_g = 2600; % Mass density of glas (kg/m3)

% Other parameters
E = 3.4e10; % Young's Modulus of Elasticity (N/m2)
Ix = 1/12*d*w^3; %(m4)
Iy = 1/12*w*d^3; %(m4)
Ib = 1/12*t*h^3; %(m4)

%Mass of column and beam
mc = w*d*rho_c; % (kg/m)
mb = h*t*rho_c; % (kg/m)
DL_story = 4E3/9.81; % (kg/m2) % ganger med E3 for at få det i N fra kN
DL_roof = 2.5E3/9.81; %(kg/m2)
LL = 0.3*2E3/9.81; % (kg/m2) % partial coefficient
t_g = 0.012; % thickness of glass (m)
t_s = 0.18; % thickness of slab [m]
m_glass = t_g*rho_g*H; % (kg/m)
m_slab = 1/2*L*W*t_s*rho_c; % (kg)

M11 = mb*(W+L)+DL_story*W*1/2*L+LL*W*1/2*L+m_glass*(W+L)+m_slab;
M44 = mb*(W+L)+DL_story*W*1/2*L+LL*W*1/2*L+m_glass*(W+L)+m_slab;
M77 = mb*(W+L)+DL_story*W*1/2*L+LL*W*1/2*L+m_glass*(W+L)+m_slab;
M1010 = mb*(W+L)+DL_roof*W*1/2*L+LL*W*1/2*L+1/2*m_glass*(W+L)+m_slab;

% M11 = mb*(W)+DL_story*W*1/2*L+m_glass*(W)+m_slab;
% M44 = mb*(W+L)+DL_story*W*1/2*L+LL*W*1/2*L+m_glass*(W+L)+m_slab;
% M77 = mb*(W+L)+DL_story*W*1/2*L+LL*W*1/2*L+m_glass*(W+L)+m_slab;
% M1010 = mb*(W)+DL_roof*W*1/2*L+1/2*m_glass*(W)+m_slab;

%% Stiffness matrix for beams and colums. 
% søjler i x
Kel_col_x = E*Ix*[12/H^3,6/H^2,-12/H^3,6/H^2;
                  6/H^2,4/H,-6/H^2,2/H;
                  -12/H^3,-6/H^2,12/H^3,-6/H^2;
                  6/H^2,2/H,-6/H^2,4/H];
              
% søjler i y
Kel_col_y = E*Iy*[12/H^3,6/H^2,-12/H^3,6/H^2;
                  6/H^2,4/H,-6/H^2,2/H;
                  -12/H^3,-6/H^2,12/H^3,-6/H^2;
                  6/H^2,2/H,-6/H^2,4/H];
    
% bjælker i x
Kel_beam_x = E*Ib*[12/W^3,6/W^2,-12/W^3,6/W^2;
                  6/W^2,4/W,-6/W^2,2/W;
                  -12/W^3,-6/W^2,12/W^3,-6/W^2;
                  6/W^2,2/W,-6/W^2,4/W];
              
% søjler i y
Kel_beam_y = E*Ib*[12/L^3,6/L^2,-12/L^3,6/L^2;
                  6/L^2,4/L,-6/L^2,2/L;
                  -12/L^3,-6/L^2,12/L^3,-6/L^2;
                  6/L^2,2/L,-6/L^2,4/L];

%% Mass matrix for beams and colums
Mel_col_x = mc*H/420*[156, 22*H, 54, -13*H;
                    22*H, 4*H^2, 13*H, -3*H^2;
                    54, 13*H, 156, -22*H;
                    -13*H,-3*H^2, -22*H, 4*H^2];

%Same as in x, just renamed.                
Mel_col_y = mc*H/420*[156, 22*H, 54, -13*H;
                    22*H, 4*H^2, 13*H, -3*H^2;
                    54, 13*H, 156, -22*H;
                    -13*H,-3*H^2, -22*H, 4*H^2];
                
Mel_beam_x = mb*W/420*[156, 22*W, 54, -13*W;
                    22*W, 4*W^2, 13*W, -3*W^2;
                    54, 13*W, 156, -22*W;
                    -13*W,-3*W^2, -22*W, 4*W^2];
                
Mel_beam_y = mb*L/420*[156, 22*L, 54, -13*L;
                    22*L, 4*L^2, 13*L, -3*L^2;
                    54, 13*L, 156, -22*L;
                    -13*L,-3*L^2, -22*L, 4*L^2];
                
%% Stiffness matix x-x
Kx = zeros(12);

Kx(1,1) = Kel_col_x(1,1)*2+Kel_col_x(3,3)*2;
Kx(2,1) = Kel_col_x(2,1)+Kel_col_x(4,3);
Kx(3,1) = Kel_col_x(2,1)+Kel_col_x(4,3);
Kx(4,1) = Kel_col_x(1,3)*2;
Kx(5,1) = Kel_col_x(2,3);
Kx(6,1) = Kel_col_x(2,3);

Kx(1,2) = Kel_col_x(1,2)+Kel_col_x(3,4);
Kx(2,2) = Kel_col_x(2,2)+Kel_col_x(4,4)+Kel_beam_x(4,4);
Kx(3,2) = Kel_beam_x(2,4); 
Kx(4,2) = Kel_col_x(1,4);
Kx(5,2) = Kel_col_x(2,4);

Kx(1,3) = Kel_col_x(1,2)+Kel_col_x(3,4);
Kx(2,3) = Kel_beam_x(4,2);
Kx(3,3) = Kel_col_x(2,2)+Kel_col_x(4,4)+Kel_beam_x(2,2);
Kx(4,3) = Kel_col_x(1,4);
Kx(6,3) = Kel_col_x(2,4);

Kx(1,4) = Kel_col_x(3,1)+Kel_col_x(3,1);
Kx(2,4) = Kel_col_x(4,1);
Kx(3,4) = Kel_col_x(4,1);
Kx(4,4) = Kel_col_x(1,1)*2+Kel_col_x(3,3)*2; 
Kx(5,4) = Kel_col_x(2,1)+Kel_col_x(4,3); 
Kx(6,4) = Kel_col_x(2,1)+Kel_col_x(4,3);
Kx(7,4) = Kel_col_x(1,3)*2; 
Kx(8,4) = Kel_col_x(2,3);
Kx(9,4) = Kel_col_x(2,3);

Kx(1,5) = Kel_col_x(3,2);
Kx(2,5) = Kel_col_x(4,2);
Kx(4,5) = Kel_col_x(1,2)+Kel_col_x(3,4);
Kx(5,5) = Kel_col_x(2,2)+Kel_col_x(4,4)+Kel_beam_x(4,4);
Kx(6,5) = Kel_beam_x(2,4);
Kx(7,5) = Kel_col_x(1,4);
Kx(8,5) = Kel_col_x(2,4);

Kx(1,6) = Kel_col_x(3,2);
Kx(3,6) = Kel_col_x(4,2);
Kx(4,6) = Kel_col_x(1,2)+Kel_col_x(3,4);
Kx(5,6) = Kel_beam_x(4,2);
Kx(6,6) = Kel_col_x(2,2)+Kel_col_x(4,4)+Kel_beam_x(2,2);
Kx(7,6) = Kel_col_x(1,4);
Kx(9,6) = Kel_col_x(2,4);

Kx(4,7) = Kel_col_x(3,1)+Kel_col_x(3,1);
Kx(5,7) = Kel_col_x(4,1);
Kx(6,7) = Kel_col_x(4,1);
Kx(7,7) = Kel_col_x(1,1)*2+Kel_col_x(3,3)*2; 
Kx(8,7) = Kel_col_x(2,1)+Kel_col_x(4,3); 
Kx(9,7) = Kel_col_x(2,1)+Kel_col_x(4,3);
Kx(10,7) = Kel_col_x(1,3)*2; 
Kx(11,7) = Kel_col_x(2,3);
Kx(12,7) = Kel_col_x(2,3);

Kx(4,8) = Kel_col_x(3,2);
Kx(5,8) = Kel_col_x(4,2);
Kx(7,8) = Kel_col_x(1,2)+Kel_col_x(3,4);
Kx(8,8) = Kel_col_x(2,2)+Kel_col_x(4,4)+Kel_beam_x(4,4);
Kx(9,8) = Kel_beam_x(2,4);
Kx(10,8) = Kel_col_x(1,4);
Kx(11,8) = Kel_col_x(2,4);

Kx(4,9) = Kel_col_x(3,2);
Kx(6,9) = Kel_col_x(4,2);
Kx(7,9) = Kel_col_x(1,2)+Kel_col_x(3,4);
Kx(8,9) = Kel_beam_x(4,2);
Kx(9,9) = Kel_col_x(2,2)+Kel_col_x(4,4)+Kel_beam_x(2,2);
Kx(10,9) = Kel_col_x(1,4);
Kx(12,9) = Kel_col_x(2,4);

Kx(7,10) = Kel_col_x(3,1)*2;
Kx(8,10) = Kel_col_x(4,1);
Kx(9,10) = Kel_col_x(4,1);
Kx(10,10) = Kel_col_x(1,1)*2;
Kx(11,10) = Kel_col_x(2,1);
Kx(12,10) = Kel_col_x(2,1);

Kx(7,11) = Kel_col_x(3,2);
Kx(8,11) = Kel_col_x(2,4);
Kx(10,11) = Kel_col_x(1,2);
Kx(11,11) = Kel_beam_x(4,4)+Kel_col_x(2,2);
Kx(12,11) = Kel_beam_x(2,4);

Kx(7,12) = Kel_col_x(3,2);
Kx(9,12) = Kel_col_x(2,4);
Kx(10,12) = Kel_col_x(1,2);
Kx(11,12) = Kel_beam_x(4,2);
Kx(12,12) = Kel_col_x(2,2)+Kel_beam_x(2,2);

Kx;

%check if it is symmetric
diag(Kx*transpose(Kx)^(-1));


%% Stiffness matrix y-y
Ky = zeros(12);

Ky(1,1) = Kel_col_y(1,1)*2+Kel_col_y(3,3)*2;
Ky(2,1) = Kel_col_y(2,1)+Kel_col_y(4,3);
Ky(3,1) = Kel_col_y(2,1)+Kel_col_y(4,3);
Ky(4,1) = Kel_col_y(1,3)*2;
Ky(5,1) = Kel_col_y(2,3);
Ky(6,1) = Kel_col_y(2,3);

Ky(1,2) = Kel_col_y(1,2)+Kel_col_y(3,4);
Ky(2,2) = Kel_col_y(2,2)+Kel_col_y(4,4)+Kel_beam_y(4,4);
Ky(3,2) = Kel_beam_y(2,4); 
Ky(4,2) = Kel_col_y(1,4);
Ky(5,2) = Kel_col_y(2,4);

Ky(1,3) = Kel_col_y(1,2)+Kel_col_y(3,4);
Ky(2,3) = Kel_beam_y(4,2);
Ky(3,3) = Kel_col_y(2,2)+Kel_col_y(4,4)+Kel_beam_y(2,2);
Ky(4,3) = Kel_col_y(1,4);
Ky(6,3) = Kel_col_y(2,4);

Ky(1,4) = Kel_col_y(3,1)+Kel_col_y(3,1);
Ky(2,4) = Kel_col_y(4,1);
Ky(3,4) = Kel_col_y(4,1);
Ky(4,4) = Kel_col_y(1,1)*2+Kel_col_y(3,3)*2; 
Ky(5,4) = Kel_col_y(2,1)+Kel_col_y(4,3); 
Ky(6,4) = Kel_col_y(2,1)+Kel_col_y(4,3);
Ky(7,4) = Kel_col_y(1,3)*2; 
Ky(8,4) = Kel_col_y(2,3);
Ky(9,4) = Kel_col_y(2,3);

Ky(1,5) = Kel_col_y(3,2);
Ky(2,5) = Kel_col_y(4,2);
Ky(4,5) = Kel_col_y(1,2)+Kel_col_y(3,4);
Ky(5,5) = Kel_col_y(2,2)+Kel_col_y(4,4)+Kel_beam_y(4,4);
Ky(6,5) = Kel_beam_y(2,4);
Ky(7,5) = Kel_col_y(1,4);
Ky(8,5) = Kel_col_y(2,4);

Ky(1,6) = Kel_col_y(3,2);
Ky(3,6) = Kel_col_y(4,2);
Ky(4,6) = Kel_col_y(1,2)+Kel_col_y(3,4);
Ky(5,6) = Kel_beam_y(4,2);
Ky(6,6) = Kel_col_y(2,2)+Kel_col_y(4,4)+Kel_beam_y(2,2);
Ky(7,6) = Kel_col_y(1,4);
Ky(9,6) = Kel_col_y(2,4);

Ky(4,7) = Kel_col_y(3,1)+Kel_col_y(3,1);
Ky(5,7) = Kel_col_y(4,1);
Ky(6,7) = Kel_col_y(4,1);
Ky(7,7) = Kel_col_y(1,1)*2+Kel_col_y(3,3)*2; 
Ky(8,7) = Kel_col_y(2,1)+Kel_col_y(4,3); 
Ky(9,7) = Kel_col_y(2,1)+Kel_col_y(4,3);
Ky(10,7) = Kel_col_y(1,3)*2; 
Ky(11,7) = Kel_col_y(2,3);
Ky(12,7) = Kel_col_y(2,3);

Ky(4,8) = Kel_col_y(3,2);
Ky(5,8) = Kel_col_y(4,2);
Ky(7,8) = Kel_col_y(1,2)+Kel_col_y(3,4);
Ky(8,8) = Kel_col_y(2,2)+Kel_col_y(4,4)+Kel_beam_y(4,4);
Ky(9,8) = Kel_beam_y(2,4);
Ky(10,8) = Kel_col_y(1,4);
Ky(11,8) = Kel_col_y(2,4);

Ky(4,9) = Kel_col_y(3,2);
Ky(6,9) = Kel_col_y(4,2);
Ky(7,9) = Kel_col_y(1,2)+Kel_col_y(3,4);
Ky(8,9) = Kel_beam_y(4,2);
Ky(9,9) = Kel_col_y(2,2)+Kel_col_y(4,4)+Kel_beam_y(2,2);
Ky(10,9) = Kel_col_y(1,4);
Ky(12,9) = Kel_col_y(2,4);

Ky(7,10) = Kel_col_y(3,1)*2;
Ky(8,10) = Kel_col_y(4,1);
Ky(9,10) = Kel_col_y(4,1);
Ky(10,10) = Kel_col_y(1,1)*2;
Ky(11,10) = Kel_col_y(2,1);
Ky(12,10) = Kel_col_y(2,1);

Ky(7,11) = Kel_col_y(3,2);
Ky(8,11) = Kel_col_y(2,4);
Ky(10,11) = Kel_col_y(1,2);
Ky(11,11) = Kel_beam_y(4,4)+Kel_col_y(2,2);
Ky(12,11) = Kel_beam_y(2,4);

Ky(7,12) = Kel_col_y(3,2);
Ky(9,12) = Kel_col_y(2,4);
Ky(10,12) = Kel_col_y(1,2);
Ky(11,12) = Kel_beam_y(4,2);
Ky(12,12) = Kel_col_y(2,2)+Kel_beam_y(2,2);

Ky; 

%check if it is symmetric
diag(Ky*transpose(Ky)^(-1));

%% Mass matrix x-x
%Missing the extra load!!!!
Mx = zeros(12);

Mx(1,1) = Mel_col_x(1,1)*2+Mel_col_x(3,3)*2+M11;
Mx(2,1) = Mel_col_x(2,1)+Mel_col_x(4,3);
Mx(3,1) = Mel_col_x(2,1)+Mel_col_x(4,3);
Mx(4,1) = Mel_col_x(1,3)*2;
Mx(5,1) = Mel_col_x(2,3);
Mx(6,1) = Mel_col_x(2,3);

Mx(1,2) = Mel_col_x(1,2)+Mel_col_x(3,4);
Mx(2,2) = Mel_col_x(2,2)+Mel_col_x(4,4)+Mel_beam_x(4,4);
Mx(3,2) = Mel_beam_x(2,4); 
Mx(4,2) = Mel_col_x(1,4);
Mx(5,2) = Mel_col_x(2,4);

Mx(1,3) = Mel_col_x(1,2)+Mel_col_x(3,4);
Mx(2,3) = Mel_beam_x(4,2);
Mx(3,3) = Mel_col_x(2,2)+Mel_col_x(4,4)+Mel_beam_x(2,2);
Mx(4,3) = Mel_col_x(1,4);
Mx(6,3) = Mel_col_x(2,4);

Mx(1,4) = Mel_col_x(3,1)+Mel_col_x(3,1);
Mx(2,4) = Mel_col_x(4,1);
Mx(3,4) = Mel_col_x(4,1);
Mx(4,4) = Mel_col_x(1,1)*2+Mel_col_x(3,3)*2+M44; 
Mx(5,4) = Mel_col_x(2,1)+Mel_col_x(4,3); 
Mx(6,4) = Mel_col_x(2,1)+Mel_col_x(4,3);
Mx(7,4) = Mel_col_x(1,3)*2; 
Mx(8,4) = Mel_col_x(2,3);
Mx(9,4) = Mel_col_x(2,3);

Mx(1,5) = Mel_col_x(3,2);
Mx(2,5) = Mel_col_x(4,2);
Mx(4,5) = Mel_col_x(1,2)+Mel_col_x(3,4);
Mx(5,5) = Mel_col_x(2,2)+Mel_col_x(4,4)+Mel_beam_x(4,4);
Mx(6,5) = Mel_beam_x(2,4);
Mx(7,5) = Mel_col_x(1,4);
Mx(8,5) = Mel_col_x(2,4);

Mx(1,6) = Mel_col_x(3,2);
Mx(3,6) = Mel_col_x(4,2);
Mx(4,6) = Mel_col_x(1,2)+Mel_col_x(3,4);
Mx(5,6) = Mel_beam_x(4,2);
Mx(6,6) = Mel_col_x(2,2)+Mel_col_x(4,4)+Mel_beam_x(2,2);
Mx(7,6) = Mel_col_x(1,4);
Mx(9,6) = Mel_col_x(2,4);

Mx(4,7) = Mel_col_x(3,1)+Mel_col_x(3,1);
Mx(5,7) = Mel_col_x(4,1);
Mx(6,7) = Mel_col_x(4,1);
Mx(7,7) = Mel_col_x(1,1)*2+Mel_col_x(3,3)*2+M77; 
Mx(8,7) = Mel_col_x(2,1)+Mel_col_x(4,3); 
Mx(9,7) = Mel_col_x(2,1)+Mel_col_x(4,3);
Mx(10,7) = Mel_col_x(1,3)*2; 
Mx(11,7) = Mel_col_x(2,3);
Mx(12,7) = Mel_col_x(2,3);

Mx(4,8) = Mel_col_x(3,2);
Mx(5,8) = Mel_col_x(4,2);
Mx(7,8) = Mel_col_x(1,2)+Mel_col_x(3,4);
Mx(8,8) = Mel_col_x(2,2)+Mel_col_x(4,4)+Mel_beam_x(4,4);
Mx(9,8) = Mel_beam_x(2,4);
Mx(10,8) = Mel_col_x(1,4);
Mx(11,8) = Mel_col_x(2,4);

Mx(4,9) = Mel_col_x(3,2);
Mx(6,9) = Mel_col_x(4,2);
Mx(7,9) = Mel_col_x(1,2)+Mel_col_x(3,4);
Mx(8,9) = Mel_beam_x(4,2);
Mx(9,9) = Mel_col_x(2,2)+Mel_col_x(4,4)+Mel_beam_x(2,2);
Mx(10,9) = Mel_col_x(1,4);
Mx(12,9) = Mel_col_x(2,4);

Mx(7,10) = Mel_col_x(3,1)*2;
Mx(8,10) = Mel_col_x(4,1);
Mx(9,10) = Mel_col_x(4,1);
Mx(10,10) = Mel_col_x(1,1)*2+M1010;
Mx(11,10) = Mel_col_x(2,1);
Mx(12,10) = Mel_col_x(2,1);

Mx(7,11) = Mel_col_x(3,2);
Mx(8,11) = Mel_col_x(2,4);
Mx(10,11) = Mel_col_x(1,2);
Mx(11,11) = Mel_beam_x(4,4)+Mel_col_x(2,2);
Mx(12,11) = Mel_beam_x(2,4);

Mx(7,12) = Mel_col_x(3,2);
Mx(9,12) = Mel_col_x(2,4);
Mx(10,12) = Mel_col_x(1,2);
Mx(11,12) = Mel_beam_x(4,2);
Mx(12,12) = Mel_col_x(2,2)+Mel_beam_x(2,2);

Mx;

%check if it is symmetric
diag(Mx*transpose(Mx)^(-1));


%% Mass matrix y-y  
%Missing the extra load!!!!
My = zeros(12);

My(1,1) = Mel_col_y(1,1)*2+Mel_col_y(3,3)*2+M11;
My(2,1) = Mel_col_y(2,1)+Mel_col_y(4,3);
My(3,1) = Mel_col_y(2,1)+Mel_col_y(4,3);
My(4,1) = Mel_col_y(1,3)*2;
My(5,1) = Mel_col_y(2,3);
My(6,1) = Mel_col_y(2,3);

My(1,2) = Mel_col_y(1,2)+Mel_col_y(3,4);
My(2,2) = Mel_col_y(2,2)+Mel_col_y(4,4)+Mel_beam_y(4,4);
My(3,2) = Mel_beam_y(2,4); 
My(4,2) = Mel_col_y(1,4);
My(5,2) = Mel_col_y(2,4);

My(1,3) = Mel_col_y(1,2)+Mel_col_y(3,4);
My(2,3) = Mel_beam_y(4,2);
My(3,3) = Mel_col_y(2,2)+Mel_col_y(4,4)+Mel_beam_y(2,2);
My(4,3) = Mel_col_y(1,4);
My(6,3) = Mel_col_y(2,4);

My(1,4) = Mel_col_y(3,1)+Mel_col_y(3,1);
My(2,4) = Mel_col_y(4,1);
My(3,4) = Mel_col_y(4,1);
My(4,4) = Mel_col_y(1,1)*2+Mel_col_y(3,3)*2+M44; 
My(5,4) = Mel_col_y(2,1)+Mel_col_y(4,3); 
My(6,4) = Mel_col_y(2,1)+Mel_col_y(4,3);
My(7,4) = Mel_col_y(1,3)*2; 
My(8,4) = Mel_col_y(2,3);
My(9,4) = Mel_col_y(2,3);

My(1,5) = Mel_col_y(3,2);
My(2,5) = Mel_col_y(4,2);
My(4,5) = Mel_col_y(1,2)+Mel_col_y(3,4);
My(5,5) = Mel_col_y(2,2)+Mel_col_y(4,4)+Mel_beam_y(4,4);
My(6,5) = Mel_beam_y(2,4);
My(7,5) = Mel_col_y(1,4);
My(8,5) = Mel_col_y(2,4);

My(1,6) = Mel_col_y(3,2);
My(3,6) = Mel_col_y(4,2);
My(4,6) = Mel_col_y(1,2)+Mel_col_y(3,4);
My(5,6) = Mel_beam_y(4,2);
My(6,6) = Mel_col_y(2,2)+Mel_col_y(4,4)+Mel_beam_y(2,2);
My(7,6) = Mel_col_y(1,4);
My(9,6) = Mel_col_y(2,4);

My(4,7) = Mel_col_y(3,1)+Mel_col_y(3,1);
My(5,7) = Mel_col_y(4,1);
My(6,7) = Mel_col_y(4,1);
My(7,7) = Mel_col_y(1,1)*2+Mel_col_y(3,3)*2+M77; 
My(8,7) = Mel_col_y(2,1)+Mel_col_y(4,3); 
My(9,7) = Mel_col_y(2,1)+Mel_col_y(4,3);
My(10,7) = Mel_col_y(1,3)*2; 
My(11,7) = Mel_col_y(2,3);
My(12,7) = Mel_col_y(2,3);

My(4,8) = Mel_col_y(3,2);
My(5,8) = Mel_col_y(4,2);
My(7,8) = Mel_col_y(1,2)+Mel_col_y(3,4);
My(8,8) = Mel_col_y(2,2)+Mel_col_y(4,4)+Mel_beam_y(4,4);
My(9,8) = Mel_beam_y(2,4);
My(10,8) = Mel_col_y(1,4);
My(11,8) = Mel_col_y(2,4);

My(4,9) = Mel_col_y(3,2);
My(6,9) = Mel_col_y(4,2);
My(7,9) = Mel_col_y(1,2)+Mel_col_y(3,4);
My(8,9) = Mel_beam_y(4,2);
My(9,9) = Mel_col_y(2,2)+Mel_col_y(4,4)+Mel_beam_y(2,2);
My(10,9) = Mel_col_y(1,4);
My(12,9) = Mel_col_y(2,4);

My(7,10) = Mel_col_y(3,1)*2;
My(8,10) = Mel_col_y(4,1);
My(9,10) = Mel_col_y(4,1);
My(10,10) = Mel_col_y(1,1)*2+M1010;
My(11,10) = Mel_col_y(2,1);
My(12,10) = Mel_col_y(2,1);

My(7,11) = Mel_col_y(3,2);
My(8,11) = Mel_col_y(2,4);
My(10,11) = Mel_col_y(1,2);
My(11,11) = Mel_beam_y(4,4)+Mel_col_y(2,2);
My(12,11) = Mel_beam_y(2,4);

My(7,12) = Mel_col_y(3,2);
My(9,12) = Mel_col_y(2,4);
My(10,12) = Mel_col_y(1,2);
My(11,12) = Mel_beam_y(4,2);
My(12,12) = Mel_col_y(2,2)+Mel_beam_y(2,2);

My;

%check if it is symmetric
diag(My*transpose(My)^(-1));

%% Frequencies for x-x 
% Mx = diag(diag(Mx))
[Ux,Dx] = eig(Kx,Mx);
omega_x = diag(sqrt(Dx)); % Natural angular frequency
Tx = 2*pi./omega_x;

%Frequencies
f_xx = 1./Tx;

%sorting 
[omega_x,isort_x] = sort(omega_x);
Ux = Ux(:,isort_x);

% Normalization
MUx = max(Ux);
mUx = min(Ux);
n_modes = 12; 
n_DOF = 12;
for j = 1:n_modes
    if abs(MUx(j)) > abs(mUx(j))
        mxUx(j) = MUx(j);
    else
        mxUx(j) = mUx(j);
    end
    for l = 1:n_DOF
        Vx(l,j) = Ux(l,j)/mxUx(j); % mode shape
    end
end

%% Frequencies for y-y
%My = diag(diag((My)))
[Uy,Dy] = eig(Ky,My);
omega_y = diag(sqrt(Dy)); % Natural angular frequency

[omega_y,isort_y] = sort(omega_y);
Uy = Uy(:,isort_y);
Ty = 2*pi./omega_y;
%Frequencies
f_yy = 1./Ty;

omega_x;
omega_y;
Tx;
Ty;

% Normalization
MUy = max(Uy);
mUy = min(Uy);
n_modes = 12; 
n_DOF = 12;
for j = 1:n_modes
    if abs(MUy(j)) > abs(mUy(j))
        myUy(j) = MUy(j);
    else
        myUy(j) = mUy(j);
    end
    for l = 1:n_DOF
        Vy(l,j) = Uy(l,j)/myUy(j); % mode shape
    end
end




% Issym(Kx)
% issym(Ky)
% issym(Mx)
% issym(My)
%% mode shape matrix
phi_x = Vx;
phi_y = Vy;
%% Modal mass [kg]
M_n_x = phi_x'*Mx*phi_x;
M_n_y = phi_y'*My*phi_y;

%% Modal stiffness [i think N/m]
K_n_x = phi_x'*Kx*phi_x;
K_n_y = phi_y'*Ky*phi_y;

%% Gamma_n [-]
r = [1,0,0,1,0,0,1,0,0,1,0,0]';
L_n_x = phi_x'*Mx*r;
L_n_y = phi_y'*My*r;

 MD_n_x = diag(M_n_x);
 MD_n_y = diag(M_n_y);
% 
% Gamma_n_x = zeros(12,1);
% Gamma_n_y = zeros(12,1);
% 
% for i = 1:12 
%     Gamma_n_x(i) = L_n_x(i)/MD_n_x(i);
%     Gamma_n_y(i) = L_n_y(i)/MD_n_y(i);
% end
Gamma_n_x = L_n_x./diag(M_n_x);
Gamma_n_y = L_n_y./diag(M_n_y);

%% Effective modal mass  [kg]
M_effn_x = L_n_x.^2./diag(M_n_x);
M_effn_y = L_n_y.^2./diag(M_n_y);


%% Find out how many modes shapes that is needed to be over the 90% threshold
m_tot=(M11+M44+M77+M1010+mc*8*H);
mass_procent_x = M_effn_x/m_tot*100; % (%)
mass_procent_y = M_effn_y/m_tot*100; % (%)
cumx = cumsum(mass_procent_x);
cumy = cumsum(mass_procent_y);
%% Topology for y-axis
% Coordinates of nodes X = [x y], 
X_x = [ 0       0
        W       0
        0       H
        W       H
        0       2*H
        W       2*H
        0       3*H
        W       3*H
        0       4*H
        W       4*H];


% Topology matrix T = [node1 node2 propno],
T_x = [   1   3   1
        2   4   1
        3   4   1
        3   5   1
        4   6   1
        5   6   1
        5   7   1
        6   8   1
        7   8   1
        7   9   1
        8   10   1
        9   10   1];

%% PLOT for x-axis
% Having the normalized mode shapes
phi_x = Vx;

modes = 12;

E1x = zeros(modes,4);
E2x = zeros(modes,4);
E3x = zeros(modes,4);
E4x = zeros(modes,4);
E5x = zeros(modes,4);
E6x = zeros(modes,4);
E7x = zeros(modes,4);
E8x = zeros(modes,4);
E9x = zeros(modes,4);
E10x = zeros(modes,4);
E11x = zeros(modes,4);
E12x = zeros(modes,4);
% Making the global DOF into the local DOF for each element

for i = 1:modes
    E1x(i,:) = [phi_x(1,i), phi_x(3,i), 0, 0];
    E2x(i,:) = [0, phi_x(3,i), 0, phi_x(2,i)];
    E3x(i,:) = [phi_x(1,i), phi_x(2,i), 0, 0];
    E4x(i,:) = [phi_x(4,i), phi_x(6,i), phi_x(1,i), phi_x(3,i)];
    E5x(i,:) = [0, phi_x(6,i), 0, phi_x(5,i)];
    E6x(i,:) = [phi_x(4,i), phi_x(5,i), phi_x(1,i), phi_x(2,i)];
    E7x(i,:) = [phi_x(7,i), phi_x(9,i), phi_x(4,i), phi_x(6,i)];
    E8x(i,:) = [0, phi_x(9,i), 0, phi_x(8,i)];
    E9x(i,:) = [phi_x(7,i), phi_x(8,i), phi_x(4,i), phi_x(5,i)];
    E10x(i,:) = [phi_x(10,i), phi_x(12,i), phi_x(7,i), phi_x(9,i)];
    E11x(i,:) = [0, phi_x(12,i), 0, phi_x(11,i)];
    E12x(i,:) = [phi_x(10,i), phi_x(11,i), phi_x(7,i), phi_x(8,i)];
end

x_col = linspace(0,H,110);

N1x = zeros(length(x_col),1);
N2x = zeros(length(x_col),1);
N3x = zeros(length(x_col),1);
N4x = zeros(length(x_col),1);
% For the columns

for i = 1:length(x_col)
    N1x(i) = 1-3*(x_col(i)/H)^2+2*(x_col(i)/H)^3;
    N2x(i) = x_col(i)*(1-x_col(i)/H)^2;
    N3x(i) = 3*(x_col(i)/H)^2-2*(x_col(i)/H)^3;
    N4x(i) = (x_col(i)/H-1)*x_col(i)^2/H;
end

x_beam = linspace(0,W,110);

N1bx = zeros(length(x_beam),1);
N2bx = zeros(length(x_beam),1);
N3bx = zeros(length(x_beam),1);
N4bx = zeros(length(x_beam),1);
% % For the beams in the x direction
for i = 1:length(x_beam)
    N1bx(i) = 1-3*(x_beam(i)/W)^2+2*(x_beam(i)/W)^3;
    N2bx(i) = x_beam(i)*(1-x_beam(i)/W)^2;
    N3bx(i) = 3*(x_beam(i)/W)^2-2*(x_beam(i)/W)^3;
    N4bx(i) = (x_beam(i)/W-1)*x_beam(i)^2/W;
end

Nx1 = zeros(modes,length(x_beam));
Nx2 = zeros(modes,length(x_beam));
Nx3 = zeros(modes,length(x_beam));
Nx4 = zeros(modes,length(x_beam));
Nx5 = zeros(modes,length(x_beam));
Nx6 = zeros(modes,length(x_beam));
Nx7 = zeros(modes,length(x_beam));
Nx8 = zeros(modes,length(x_beam));
Nx9 = zeros(modes,length(x_beam));
Nx10 = zeros(modes,length(x_beam));
Nx11 = zeros(modes,length(x_beam));
Nx12 = zeros(modes,length(x_beam));

for i = 1:modes
    for j = 1:length(x_col)
        Nx1(i,:) = (E1x(i,1)*N1x'+E1x(i,2)*N2x'+E1x(i,3)*N3x'+E1x(i,4)*N4x');
        Nx2(i,:) = (E2x(i,1)*N1bx'+E2x(i,2)*N2bx'+E2x(i,3)*N3bx'+E2x(i,4)*N4bx');
        Nx3(i,:) = (E3x(i,1)*N1x'+E3x(i,2)*N2x'+E3x(i,3)*N3x'+E3x(i,4)*N4x');
        Nx4(i,:) = (E4x(i,1)*N1x'+E4x(i,2)*N2x'+E4x(i,3)*N3x'+E4x(i,4)*N4x');
        Nx5(i,:) = (E5x(i,1)*N1bx'+E5x(i,2)*N2bx'+E5x(i,3)*N3bx'+E5x(i,4)*N4bx');
        Nx6(i,:) = (E6x(i,1)*N1x'+E6x(i,2)*N2x'+E6x(i,3)*N3x'+E6x(i,4)*N4x');
        Nx7(i,:) = (E7x(i,1)*N1x'+E7x(i,2)*N2x'+E7x(i,3)*N3x'+E7x(i,4)*N4x');
        Nx8(i,:) = (E8x(i,1)*N1bx'+E8x(i,2)*N2bx'+E8x(i,3)*N3bx'+E8x(i,4)*N4bx');
        Nx9(i,:) = (E9x(i,1)*N1x'+E9x(i,2)*N2x'+E9x(i,3)*N3x'+E9x(i,4)*N4x');
        Nx10(i,:) = (E10x(i,1)*N1x'+E10x(i,2)*N2x'+E10x(i,3)*N3x'+E10x(i,4)*N4x');
        Nx11(i,:) = (E11x(i,1)*N1bx'+E11x(i,2)*N2bx'+E11x(i,3)*N3bx'+E11x(i,4)*N4bx');
        Nx12(i,:) = (E12x(i,1)*N1x'+E12x(i,2)*N2x'+E12x(i,3)*N3x'+E12x(i,4)*N4x');
    end
end

k=zeros(4,1);
k2=zeros(4,1);
k3=zeros(4,1);
h=zeros(4,1);
h2=zeros(4,1);
h3=zeros(4,1);
fig = figure;
for i = 1:modes
    for j = 1:4
% rotation of 270 degrees
        h(j) = 0; % Translation in x
        k(j) = j*H; % Translation in y
    end

        a = 3/2 * pi;
        A1 = [cos(a) -sin(a) h(1);
                sin(a) cos(a)  k(1);
                0      0       1];
        A2 = [cos(a) -sin(a) h(2);
                sin(a) cos(a)  k(2);
                0      0       1];
        A3 = [cos(a) -sin(a) h(3);
                sin(a) cos(a)  k(3);
                0      0       1];
        A4 = [cos(a) -sin(a) h(4);
                sin(a) cos(a)  k(4);
                0      0       1];

        

        beam1 = [x_col;Nx1(i,:);ones(1,length(x_beam))];
        beam4 = [x_col;Nx4(i,:);ones(1,length(x_beam))];
        beam7 = [x_col;Nx7(i,:);ones(1,length(x_beam))];
        beam10 = [x_col;Nx10(i,:);ones(1,length(x_beam))];
    
    
        % Tranforming the beam
        beam_tranf1 = A1 * beam1;
        beam_tranf4 = A2 * beam4;
        beam_tranf7 = A3 * beam7;
        beam_tranf10 = A4 * beam10;

    for q = 1:4
        % rotation of 270 degrees
        h2(q) = 0+phi_x(1+(q-1)*3,i); % Translation in x ændrer 1 til i nå der plottes
        k2(q) = q*H; % Translation in y
    end

        a = 0 * pi;
        A5 = [cos(a) -sin(a) h2(1);
             sin(a) cos(a)  k2(1);
             0      0       1];
        A6 = [cos(a) -sin(a) h2(2);
             sin(a) cos(a)  k2(2);
             0      0       1];      
        A7 = [cos(a) -sin(a) h2(3);
             sin(a) cos(a)  k2(3);
             0      0       1];
        A8 = [cos(a) -sin(a) h2(4);
             sin(a) cos(a)  k2(4);
             0      0       1];
 
         beam2 = [x_beam;Nx2(i,:);ones(1,length(x_beam))];
         beam5 = [x_beam;Nx5(i,:);ones(1,length(x_beam))];
         beam8 = [x_beam;Nx8(i,:);ones(1,length(x_beam))];
         beam11 = [x_beam;Nx11(i,:);ones(1,length(x_beam))];
         
         beam_tranf2 = A5 * beam2;
         beam_tranf5 = A6 * beam5;
         beam_tranf8 = A7 * beam8;
         beam_tranf11 = A8 * beam11;
  
    for n = 1:4
        % rotation of 270 degrees
        h3(n) = W; % Translation in x
        k3(n) = n*H; % Translation in y
    end
        a = 3/2 * pi;
        A9 = [cos(a) -sin(a) h3(1);
             sin(a) cos(a)  k3(1);
             0      0       1];
        A10 = [cos(a) -sin(a) h3(2);
             sin(a) cos(a)  k3(2);
             0      0       1];
        A11 = [cos(a) -sin(a) h3(3);
             sin(a) cos(a)  k3(3);
             0      0       1];
        A12 = [cos(a) -sin(a) h3(4);
             sin(a) cos(a)  k3(4);
             0      0       1];
    
         
       beam3 = [x_col;Nx3(i,:);ones(1,length(x_beam))];
       beam6 = [x_col;Nx6(i,:);ones(1,length(x_beam))];
       beam9 = [x_col;Nx9(i,:);ones(1,length(x_beam))];
       beam12 = [x_col;Nx12(i,:);ones(1,length(x_beam))];
       
       beam_tranf3 = A9 * beam3;
       beam_tranf6 = A10 * beam6;
       beam_tranf9 = A11 * beam9;
       beam_tranf12 = A12 * beam12;
       
% figure ()
% subplot(3,4,i)
% hold on 
% plot(beam_tranf1(1,:),beam_tranf1(2,:))
% plot(beam_tranf2(1,:),beam_tranf2(2,:))
% plot(beam_tranf3(1,:),beam_tranf3(2,:))
% plot(beam_tranf4(1,:),beam_tranf4(2,:))
% plot(beam_tranf5(1,:),beam_tranf5(2,:))
% plot(beam_tranf6(1,:),beam_tranf6(2,:))
% plot(beam_tranf7(1,:),beam_tranf7(2,:))
% plot(beam_tranf8(1,:),beam_tranf8(2,:))
% plot(beam_tranf9(1,:),beam_tranf9(2,:))
% plot(beam_tranf10(1,:),beam_tranf10(2,:))
% plot(beam_tranf11(1,:),beam_tranf11(2,:))
% plot(beam_tranf12(1,:),beam_tranf12(2,:))
% hold off
% xlabel('Deflection ')
% ylabel('Height [m]')
% title(sprintf('Mode shape %d',i))




subplot(3,4,i)
hold on 
plot(beam_tranf1(1,:),beam_tranf1(2,:))
plot(beam_tranf2(1,:),beam_tranf2(2,:))
plot(beam_tranf3(1,:),beam_tranf3(2,:))
plot(beam_tranf4(1,:),beam_tranf4(2,:))
plot(beam_tranf5(1,:),beam_tranf5(2,:))
plot(beam_tranf6(1,:),beam_tranf6(2,:))
plot(beam_tranf7(1,:),beam_tranf7(2,:))
plot(beam_tranf8(1,:),beam_tranf8(2,:))
plot(beam_tranf9(1,:),beam_tranf9(2,:))
plot(beam_tranf10(1,:),beam_tranf10(2,:))
plot(beam_tranf11(1,:),beam_tranf11(2,:))
plot(beam_tranf12(1,:),beam_tranf12(2,:))
plotelem(T_x(:,1:3),X_x,'k:');
ylim([0 4*H+H/2])
hold off
title(sprintf('Mode shape %d',i))

end


han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]');
xlabel(han,'Deflection [-]');

%saveas(gcf,'C:\Users\Frede\Danmarks Tekniske Universitet\Katarina Kolding Skaarup - 41967 Seismic and wind engineering\Figures\LR_modeshapes_x.eps','epsc2')

%% Topology y-axis
% Coordinates of nodes X = [x y], 
X_y = [   0       0
        L       0
        0       H
        L       H
        0       2*H
        L       2*H
        0       3*H
        L       3*H
        0       4*H
        L       4*H ];


% Topology matrix T = [node1 node2 propno],
T_y = [   1   3   1
        2   4   1
        3   4   1
        3   5   1
        4   6   1
        5   6   1
        5   7   1
        6   8   1
        7   8   1
        7   9   1
        8   10   1
        9   10   1];

%% PLot for y axis
phi_y = Vy;

modes = 12;

E1y = zeros(modes,4);
E2y = zeros(modes,4);
E3y = zeros(modes,4);
E4y = zeros(modes,4);
E5y = zeros(modes,4);
E6y = zeros(modes,4);
E7y = zeros(modes,4);
E8y = zeros(modes,4);
E9y = zeros(modes,4);
E10y = zeros(modes,4);
E11y = zeros(modes,4);
E12y = zeros(modes,4);
% Making the global DOF into the local DOF for each element

for i = 1:modes
    E1y(i,:) = [phi_y(1,i), phi_y(3,i), 0, 0];
    E2y(i,:) = [0, phi_y(3,i), 0, phi_y(2,i)];
    E3y(i,:) = [phi_y(1,i), phi_y(2,i), 0, 0];
    E4y(i,:) = [phi_y(4,i), phi_y(6,i), phi_y(1,i), phi_y(3,i)];
    E5y(i,:) = [0, phi_y(6,i), 0, phi_y(5,i)];
    E6y(i,:) = [phi_y(4,i), phi_y(5,i), phi_y(1,i), phi_y(2,i)];
    E7y(i,:) = [phi_y(7,i), phi_y(9,i), phi_y(4,i), phi_y(6,i)];
    E8y(i,:) = [0, phi_y(9,i), 0, phi_y(8,i)];
    E9y(i,:) = [phi_y(7,i), phi_y(8,i), phi_y(4,i), phi_y(5,i)];
    E10y(i,:) = [phi_y(10,i), phi_y(12,i), phi_y(7,i), phi_y(9,i)];
    E11y(i,:) = [0, phi_y(12,i), 0, phi_y(11,i)];
    E12y(i,:) = [phi_y(10,i), phi_y(11,i), phi_y(7,i), phi_y(8,i)];
end

y_col = linspace(0,H,110);

N1y = zeros(length(y_col),1);
N2y = zeros(length(y_col),1);
N3y = zeros(length(y_col),1);
N4y = zeros(length(y_col),1);
% For the columns

for i = 1:length(y_col)
    N1y(i) = 1-3*(y_col(i)/H)^2+2*(y_col(i)/H)^3;
    N2y(i) = y_col(i)*(1-y_col(i)/H)^2;
    N3y(i) = 3*(y_col(i)/H)^2-2*(y_col(i)/H)^3;
    N4y(i) = (y_col(i)/H-1)*y_col(i)^2/H;
end

y_beam = linspace(0,L,110);

N1by = zeros(length(y_beam),1);
N2by = zeros(length(y_beam),1);
N3by = zeros(length(y_beam),1);
N4by = zeros(length(y_beam),1);
% % For the beams in the y direction
for i = 1:length(y_beam)
    N1by(i) = 1-3*(y_beam(i)/L)^2+2*(y_beam(i)/L)^3;
    N2by(i) = y_beam(i)*(1-y_beam(i)/L)^2;
    N3by(i) = 3*(y_beam(i)/L)^2-2*(y_beam(i)/L)^3;
    N4by(i) = (y_beam(i)/L-1)*y_beam(i)^2/L;
end

Ny1 = zeros(modes,length(y_beam));
Ny2 = zeros(modes,length(y_beam));
Ny3 = zeros(modes,length(y_beam));
Ny4 = zeros(modes,length(y_beam));
Ny5 = zeros(modes,length(y_beam));
Ny6 = zeros(modes,length(y_beam));
Ny7 = zeros(modes,length(y_beam));
Ny8 = zeros(modes,length(y_beam));
Ny9 = zeros(modes,length(y_beam));
Ny10 = zeros(modes,length(y_beam));
Ny11 = zeros(modes,length(y_beam));
Ny12 = zeros(modes,length(y_beam));

for i = 1:modes
    for j = 1:length(y_col)
        Ny1(i,:) = (E1y(i,1)*N1y'+E1y(i,2)*N2y'+E1y(i,3)*N3y'+E1y(i,4)*N4y');
        Ny2(i,:) = (E2y(i,1)*N1by'+E2y(i,2)*N2by'+E2y(i,3)*N3by'+E2y(i,4)*N4by');
        Ny3(i,:) = (E3y(i,1)*N1y'+E3y(i,2)*N2y'+E3y(i,3)*N3y'+E3y(i,4)*N4y');
        Ny4(i,:) = (E4y(i,1)*N1y'+E4y(i,2)*N2y'+E4y(i,3)*N3y'+E4y(i,4)*N4y');
        Ny5(i,:) = (E5y(i,1)*N1by'+E5y(i,2)*N2by'+E5y(i,3)*N3by'+E5y(i,4)*N4by');
        Ny6(i,:) = (E6y(i,1)*N1y'+E6y(i,2)*N2y'+E6y(i,3)*N3y'+E6y(i,4)*N4y');
        Ny7(i,:) = (E7y(i,1)*N1y'+E7y(i,2)*N2y'+E7y(i,3)*N3y'+E7y(i,4)*N4y');
        Ny8(i,:) = (E8y(i,1)*N1by'+E8y(i,2)*N2by'+E8y(i,3)*N3by'+E8y(i,4)*N4by');
        Ny9(i,:) = (E9y(i,1)*N1y'+E9y(i,2)*N2y'+E9y(i,3)*N3y'+E9y(i,4)*N4y');
        Ny10(i,:) = (E10y(i,1)*N1y'+E10y(i,2)*N2y'+E10y(i,3)*N3y'+E10y(i,4)*N4y');
        Ny11(i,:) = (E11y(i,1)*N1by'+E11y(i,2)*N2by'+E11y(i,3)*N3by'+E11y(i,4)*N4by');
        Ny12(i,:) = (E12y(i,1)*N1y'+E12y(i,2)*N2y'+E12y(i,3)*N3y'+E12y(i,4)*N4y');
    end
end

k=zeros(4,1);
k2=zeros(4,1);
k3=zeros(4,1);
h=zeros(4,1);
h2=zeros(4,1);
h3=zeros(4,1);
fig = figure;
for i = 1:modes
    for j = 1:4
% rotation of 270 degrees
        h(j) = 0; % Translation in x
        k(j) = j*H; % Translation in y
    end

        a = 3/2 * pi;
        A1 = [cos(a) -sin(a) h(1);
                sin(a) cos(a)  k(1);
                0      0       1];
        A2 = [cos(a) -sin(a) h(2);
                sin(a) cos(a)  k(2);
                0      0       1];
        A3 = [cos(a) -sin(a) h(3);
                sin(a) cos(a)  k(3);
                0      0       1];
        A4 = [cos(a) -sin(a) h(4);
                sin(a) cos(a)  k(4);
                0      0       1];

        

        beam1 = [y_col;Ny1(i,:);ones(1,length(y_beam))];
        beam4 = [y_col;Ny4(i,:);ones(1,length(y_beam))];
        beam7 = [y_col;Ny7(i,:);ones(1,length(y_beam))];
        beam10 = [y_col;Ny10(i,:);ones(1,length(y_beam))];
    
    
        % Tranforming the beam
        beam_tranf1 = A1 * beam1;
        beam_tranf4 = A2 * beam4;
        beam_tranf7 = A3 * beam7;
        beam_tranf10 = A4 * beam10;

    for q = 1:4
        % rotation of 270 degrees
        h2(q) = 0+phi_y(1+(q-1)*3,i); % Translation in x ændrer 1 til i nå der plottes
        k2(q) = q*H; % Translation in y
    end

        a = 0 * pi;
        A5 = [cos(a) -sin(a) h2(1);
             sin(a) cos(a)  k2(1);
             0      0       1];
        A6 = [cos(a) -sin(a) h2(2);
             sin(a) cos(a)  k2(2);
             0      0       1];      
        A7 = [cos(a) -sin(a) h2(3);
             sin(a) cos(a)  k2(3);
             0      0       1];
        A8 = [cos(a) -sin(a) h2(4);
             sin(a) cos(a)  k2(4);
             0      0       1];
 
         beam2 = [y_beam;Ny2(i,:);ones(1,length(y_beam))];
         beam5 = [y_beam;Ny5(i,:);ones(1,length(y_beam))];
         beam8 = [y_beam;Ny8(i,:);ones(1,length(y_beam))];
         beam11 = [y_beam;Ny11(i,:);ones(1,length(y_beam))];
         
         beam_tranf2 = A5 * beam2;
         beam_tranf5 = A6 * beam5;
         beam_tranf8 = A7 * beam8;
         beam_tranf11 = A8 * beam11;
  
    for n = 1:4
        % rotation of 270 degrees
        h3(n) = L; % Translation in x
        k3(n) = n*H; % Translation in y
    end
        a = 3/2 * pi;
        A9 = [cos(a) -sin(a) h3(1);
             sin(a) cos(a)  k3(1);
             0      0       1];
        A10 = [cos(a) -sin(a) h3(2);
             sin(a) cos(a)  k3(2);
             0      0       1];
        A11 = [cos(a) -sin(a) h3(3);
             sin(a) cos(a)  k3(3);
             0      0       1];
        A12 = [cos(a) -sin(a) h3(4);
             sin(a) cos(a)  k3(4);
             0      0       1];
    
         
       beam3 = [y_col;Ny3(i,:);ones(1,length(y_beam))];
       beam6 = [y_col;Ny6(i,:);ones(1,length(y_beam))];
       beam9 = [y_col;Ny9(i,:);ones(1,length(y_beam))];
       beam12 = [y_col;Ny12(i,:);ones(1,length(y_beam))];
       
       beam_tranf3 = A9 * beam3;
       beam_tranf6 = A10 * beam6;
       beam_tranf9 = A11 * beam9;
       beam_tranf12 = A12 * beam12;
       
% figure ()
% subplot(3,4,i)
% hold on 
% plot(beam_tranf1(1,:),beam_tranf1(2,:))
% plot(beam_tranf2(1,:),beam_tranf2(2,:))
% plot(beam_tranf3(1,:),beam_tranf3(2,:))
% plot(beam_tranf4(1,:),beam_tranf4(2,:))
% plot(beam_tranf5(1,:),beam_tranf5(2,:))
% plot(beam_tranf6(1,:),beam_tranf6(2,:))
% plot(beam_tranf7(1,:),beam_tranf7(2,:))
% plot(beam_tranf8(1,:),beam_tranf8(2,:))
% plot(beam_tranf9(1,:),beam_tranf9(2,:))
% plot(beam_tranf10(1,:),beam_tranf10(2,:))
% plot(beam_tranf11(1,:),beam_tranf11(2,:))
% plot(beam_tranf12(1,:),beam_tranf12(2,:))
% hold off
% xlabel('Deflection ')
% ylabel('Height [m]')
% title(sprintf('Mode shape %d',i))




subplot(3,4,i)
hold on 
plot(beam_tranf1(1,:),beam_tranf1(2,:))
plot(beam_tranf2(1,:),beam_tranf2(2,:))
plot(beam_tranf3(1,:),beam_tranf3(2,:))
plot(beam_tranf4(1,:),beam_tranf4(2,:))
plot(beam_tranf5(1,:),beam_tranf5(2,:))
plot(beam_tranf6(1,:),beam_tranf6(2,:))
plot(beam_tranf7(1,:),beam_tranf7(2,:))
plot(beam_tranf8(1,:),beam_tranf8(2,:))
plot(beam_tranf9(1,:),beam_tranf9(2,:))
plot(beam_tranf10(1,:),beam_tranf10(2,:))
plot(beam_tranf11(1,:),beam_tranf11(2,:))
plot(beam_tranf12(1,:),beam_tranf12(2,:))
plotelem(T_y(:,1:3),X_y,'k:');
ylim([0 4*H+H/2])
hold off
title(sprintf('Mode shape %d',i))

end


han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]');
xlabel(han,'Deflection [-]');
%saveas(gcf,'C:\Users\Frede\Danmarks Tekniske Universitet\Katarina Kolding Skaarup - 41967 Seismic and wind engineering\Figures\LR_modeshapes_y.eps','epsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finding the mass that is used:

%% PLOT

% Coordinates of nodes X = [x y], 
X = [   0       0
        W       0
        0       H
        W       H
        0       2*H
        W       2*H
        0       3*H
        W       3*H
        0       4*H
        W       4*H ];


% Topology matrix T = [node1 node2 propno],
T = [   1   3   1
        2   4   1
        3   4   1
        3   5   1
        4   6   1
        5   6   1
        5   7   1
        6   8   1
        7   8   1
        7   9   1
        8   10   1
        9   10   1];

% Plot initial geometry
figure ()
%clf
plotelem(T(:,1:3),X,'b--');
%axis(plotaxes);
%title('Initial geometry');


%% Exercise 2 Modal Response spectrum analysis
q = 4; 

% For soil type D
S = 1.35;
TB = 0.2;
TC = 0.8; 
TD = 2.0; 
ag = 1.2*0.35*9.81; % m/s^2 remember partial coefficient
xi = 0.03; % Damping ratio in %
eta = sqrt(10/(5+100*xi));

T = linspace(0,4,1000);
Se = zeros(1,1000);

% Elastic spectrum Resonse
for i=1:1000
    if T(i) <= TB 
        Se(i) = ag*S*(1+(T(i)/TB)*(eta*2.5-1));
    elseif TB < T(i) && T(i) <= TC
        Se(i) = ag*S*eta*2.5;
    elseif TC < T(i) && T(i) <= TD
        Se(i) = ag*S*eta*2.5*(TC/T(i));
    else
        Se(i) = ag*S*eta*2.5*((TC*TD)/T(i)^2);
    end
end

%Design spectrum response
Sd = zeros(1,1000);
for i=1:1000
    if T(i) <= TB 
        Sd(i) = ag*S*(2/3+(T(i)/TB)*((2.5/q)-(2/3)));
    elseif TB < T(i) && T(i) <= TC
        Sd(i) = ag*S*(2.5/q);
    elseif TC < T(i) && T(i) <= TD
        Sd(i) = ag*S*(2.5/q)*(TC/T(i));
    else
        Sd(i) = ag*S*(2.5/q)*((TC*TD)/T(i)^2);
    end
end
fig=figure()
fig.Position=[100 100 900 550]
set(gcf,'renderer','Painters')
plot(T,Se,'-b')
hold on
plot(T,Sd,'--m')
legend('Elastic','Design','FontSize', 14)
xlabel('Period, T [s]','FontSize', 14)
ylabel('Spectral Acceleration S_e [m/s^2]','FontSize', 14)
saveas(gcf,'C:\Users\Frede\Danmarks Tekniske Universitet\Katarina Kolding Skaarup - 41967 Seismic and wind engineering\Figures\LR_spectral.eps','epsc')

% finding the spectral accelration:
bbb = zeros(12,4);
for i=1:12
bbb(i,1) = Se(find(T>Tx(i),1));
bbb(i,2) = Se(find(T>Ty(i),1));
bbb(i,3) = Sd(find(T>Tx(i),1));
bbb(i,4) = Sd(find(T>Ty(i),1));
end
%% Finding the spectral displacement in x_axis
modes = 12;
% Elastic
S_De_x = zeros(modes,1);
for i = 1:modes
S_De_x(i) = Se(find(T>Tx(i),1))*(Tx(i)/(2*pi))^2;
end

% Design
S_Dd_x = zeros(modes,1);
for i = 1:modes
S_Dd_x(i) = Sd(find(T>Tx(i),1))*(Tx(i)/(2*pi))^2;
end

%% Finding the spectral displacement in y_axis
modes = 12;
% Elastic
S_De_y = zeros(modes,1);
for i = 1:modes
S_De_y(i) = Se(find(T>Ty(i),1))*(Ty(i)/(2*pi))^2;
end

% Design
S_Dd_y = zeros(modes,1);
for i = 1:modes
S_Dd_y(i) = Sd(find(T>Ty(i),1))*(Ty(i)/(2*pi))^2;
end


%% Finding the modal displacement for x-axis
% number of modes:
modes = 12;

% The elstic spectrum [mm]
u_e_x = zeros(modes,modes);
for i = 1:modes
    u_e_x(:,i) = phi_x(:,i)*Gamma_n_x(i)*S_De_x(i);
end

% Design spectrum
u_d_x = zeros(modes,1);
for i = 1:modes
    u_d_x(:,i) = phi_x(:,i)*Gamma_n_x(i)*S_Dd_x(i);
end

%% Finding the modal displacement for y-axis
% number of modes:
modes = 12;

% The elstic spectrum
u_e_y = zeros(modes,modes);
for i = 1:modes
    u_e_y(:,i) = phi_y(:,i)*Gamma_n_y(i)*S_De_y(i);
end

% Design spectrum
u_d_y = zeros(modes,modes);
for i = 1:modes
    u_d_y(:,i) = phi_y(:,i)*Gamma_n_y(i)*S_Dd_y(i);
end

%% Plot the displacement

%% Caclulating the total deformation with the SRSS method elastic
% x-axis 
DOF = 12;
modes = 12;
r_tot_xe1 = zeros(DOF,1);
r_tot_xe2 = zeros(modes,1);
for i = 1:DOF
    for j = 1:modes
        r_tot_xe2(j) = u_e_x(i,j)^2;
        r_tot_xe1(i) = r_tot_xe1(i)+r_tot_xe2(j);
    end
end
r_tot_xe = sqrt(r_tot_xe1);

% y-axis
DOF = 12;
modes = 12;
r_tot_ye1 = zeros(DOF,1);
r_tot_ye2 = zeros(modes,1);
for i = 1:DOF
    for j = 1:modes
        r_tot_ye2(j) = u_e_y(i,j)^2;
        r_tot_ye1(i) = r_tot_ye1(i)+r_tot_ye2(j);
    end
end
r_tot_ye = sqrt(r_tot_ye1);

%% Caclulating the total deformation with the SRSS method design
% x-axis 
DOF = 12;
modes = 12;
r_tot_xd1 = zeros(DOF,1);
r_tot_xd2 = zeros(modes,1);
for i = 1:DOF
    for j = 1:modes
        r_tot_xd2(j) = u_d_x(i,j)^2;
        r_tot_xd1(i) = r_tot_xd1(i)+r_tot_xd2(j);
    end
end
r_tot_xd = sqrt(r_tot_xd1);

% y-axis
DOF = 12;
modes = 12;
r_tot_yd1 = zeros(DOF,1);
r_tot_yd2 = zeros(modes,1);
for i = 1:DOF
    for j = 1:modes
        r_tot_yd2(j) = u_d_y(i,j)^2;
        r_tot_yd1(i) = r_tot_yd1(i)+r_tot_yd2(j);
    end
end
r_tot_yd = sqrt(r_tot_yd1);

%% Caclulating the total deformation with the CQC method elstic
DOF = 12;
modes = 12;
Beta_x = zeros(DOF,modes);
Beta_y = zeros(DOF,modes);
rho_x = zeros(DOF,modes);
rho_y = zeros(DOF,modes);
for i = 1:DOF
    for j=1:modes
        Beta_x(i,j)=omega_x(i)/omega_x(j);
        Beta_y(i,j)=omega_y(i)/omega_y(j);
        rho_x(i,j)=(8*(xi^2)*(1+Beta_x(i,j)*(Beta_x(i,j)^(3/2))))/((1-Beta_x(i,j)^2)^2+4*xi^2*Beta_x(i,j)*(1+Beta_x(i,j))^2);
        rho_y(i,j)=(8*(xi^2)*(1+Beta_y(i,j)*Beta_y(i,j)^(3/2)))/((1-Beta_y(i,j)^2)^2+4*xi^2*Beta_y(i,j)*(1+Beta_y(i,j))^2);
    end
end

% x-axis 
DOF = 12;
modes = 12;
r_tot_xeCQ1 = zeros(DOF,1);

for i = 1:DOF
    for j = 1:modes
        for n = 1:modes
            if n ~= j
                r_tot_xeCQ1(i) = r_tot_xeCQ1(i)+rho_x(j,n)*u_e_x(i,j)*u_e_x(i,n);
            else 
                r_tot_xeCQ1(i) = r_tot_xeCQ1(i);
            end
        end
    end
end


r_tot_xeCQ = sqrt(r_tot_xe1+r_tot_xeCQ1);


% For y-axis
r_tot_yeCQ1 = zeros(DOF,1);
for i = 1:DOF
    for j = 1:modes
        for n = 1:modes
            if n ~= j
                r_tot_yeCQ1(i) = r_tot_yeCQ1(i)+rho_y(j,n)*u_e_y(i,j)*u_e_y(i,n);
            else 
                r_tot_yeCQ1(i) = r_tot_yeCQ1(i);
            end
        end
    end
end
r_tot_yeCQ = sqrt(r_tot_ye1+r_tot_yeCQ1);
r_tot_sam1 = [r_tot_xe,r_tot_xeCQ,r_tot_ye,r_tot_yeCQ]
%% Caclulating the total deformation with the CQC method design
% x-axis 
DOF = 12;
modes = 12;
r_tot_xdCQ1 = zeros(DOF,1);
for i = 1:DOF
    for j = 1:modes
        for n = 1:modes
            if n ~= j
                r_tot_xdCQ1(i) = r_tot_xdCQ1(i)+rho_x(j,n)*u_d_x(i,j)*u_d_x(i,n);
            else 
                r_tot_xdCQ1(i) = r_tot_xdCQ1(i);
            end
        end
    end
end


r_tot_xdCQ = sqrt(r_tot_xd1+r_tot_xdCQ1);

% For y-axis
r_tot_ydCQ1 = zeros(DOF,1);
for i = 1:DOF
    for j = 1:modes
        for n = 1:modes
            if n ~= j
                r_tot_ydCQ1(i) = r_tot_ydCQ1(i)+rho_y(j,n)*u_d_y(i,j)*u_d_y(i,n);
            else 
                r_tot_ydCQ1(i) = r_tot_ydCQ1(i);
            end
        end
    end
end
r_tot_ydCQ = sqrt(r_tot_yd1+r_tot_ydCQ1);
r_tot_sam2 = [r_tot_xd,r_tot_xdCQ,r_tot_yd,r_tot_ydCQ];
totall = [r_tot_xe,r_tot_xeCQ,r_tot_xd,r_tot_xdCQ,r_tot_ye,r_tot_yeCQ,r_tot_yd,r_tot_ydCQ]*1000;

%% Calculating the sectional forces 
% x-axis

% The dispclacement related to the design spectrum:
Uxx1 = zeros(4,modes);
Uxx2 = zeros(4,modes);
Uxx3 = zeros(4,modes);
Uxx4 = zeros(4,modes);
Uxx5 = zeros(4,modes);
Uxx6 = zeros(4,modes);
Uxx7 = zeros(4,modes);
Uxx8 = zeros(4,modes);
Uxx9 = zeros(4,modes);
Uxx10 = zeros(4,modes);
Uxx11 = zeros(4,modes);
Uxx12 = zeros(4,modes);

% oppe fra og ned
for i = 1:modes
    Uxx1(:,i) = [u_d_x(1,i), u_d_x(3,i), 0, 0];
    Uxx2(:,i) = [0, u_d_x(3,i), 0, u_d_x(2,i)];
    Uxx3(:,i) = [u_d_x(1,i), u_d_x(2,i), 0, 0];
    Uxx4(:,i) = [u_d_x(4,i), u_d_x(6,i), u_d_x(1,i), u_d_x(3,i)];
    Uxx5(:,i) = [0, u_d_x(6,i), 0, u_d_x(5,i)];
    Uxx6(:,i) = [u_d_x(4,i), u_d_x(5,i), u_d_x(1,i), u_d_x(2,i)];
    Uxx7(:,i) = [u_d_x(7,i), u_d_x(9,i), u_d_x(4,i), u_d_x(6,i)];
    Uxx8(:,i) = [0, u_d_x(9,i), 0, u_d_x(8,i)];
    Uxx9(:,i) = [u_d_x(7,i), u_d_x(8,i), u_d_x(4,i), u_d_x(5,i)];
    Uxx10(:,i) = [u_d_x(10,i), u_d_x(12,i), u_d_x(7,i), u_d_x(9,i)];
    Uxx11(:,i) = [0, u_d_x(12,i), 0, u_d_x(11,i)];
    Uxx12(:,i) = [u_d_x(10,i), u_d_x(11,i), u_d_x(7,i), u_d_x(8,i)];
end
% ned fra og op
% for i = 1:modes
%     Uxx1(:,i) = [0, 0, u_d_x(1,i), u_d_x(3,i)];
%     Uxx2(:,i) = [0, u_d_x(3,i), 0, u_d_x(2,i)];
%     Uxx3(:,i) = [0, 0, u_d_x(1,i), u_d_x(2,i)];
%     Uxx4(:,i) = [u_d_x(1,i), u_d_x(3,i), u_d_x(4,i), u_d_x(6,i)];
%     Uxx5(:,i) = [0, u_d_x(6,i), 0, u_d_x(5,i)];
%     Uxx6(:,i) = [u_d_x(1,i), u_d_x(2,i), u_d_x(4,i), u_d_x(5,i)];
%     Uxx7(:,i) = [u_d_x(4,i), u_d_x(6,i), u_d_x(7,i), u_d_x(9,i)];
%     Uxx8(:,i) = [0, u_d_x(9,i), 0, u_d_x(8,i)];
%     Uxx9(:,i) = [u_d_x(4,i), u_d_x(5,i), u_d_x(7,i), u_d_x(8,i)];
%     Uxx10(:,i) = [u_d_x(7,i), u_d_x(9,i), u_d_x(10,i), u_d_x(12,i)];
%     Uxx11(:,i) = [0, u_d_x(12,i), 0, u_d_x(11,i)];
%     Uxx12(:,i) = [u_d_x(7,i), u_d_x(8,i), u_d_x(10,i), u_d_x(11,i)];
% end

% The local stiffness matrix for coloumn
K_column_x = E*Ix*[12/H^3,6/H^2,-12/H^3,6/H^2;
                  6/H^2,4/H,-6/H^2,2/H;
                  -12/H^3,-6/H^2,12/H^3,-6/H^2;
                  6/H^2,2/H,-6/H^2,4/H];

% The local stiffness matrix for beam
K_beam_x = E*Ib*[12/W^3,6/W^2,-12/W^3,6/W^2;
                  6/W^2,4/W,-6/W^2,2/W;
                  -12/W^3,-6/W^2,12/W^3,-6/W^2;
                  6/W^2,2/W,-6/W^2,4/W];


% Sectional forces for each element is found:
FxE1 = zeros(4,modes);
FxE2 = zeros(4,modes);
FxE3 = zeros(4,modes);
FxE4 = zeros(4,modes);
FxE5 = zeros(4,modes);
FxE6 = zeros(4,modes);
FxE7 = zeros(4,modes);
FxE8 = zeros(4,modes);
FxE9 = zeros(4,modes);
FxE10 = zeros(4,modes);
FxE11 = zeros(4,modes);
FxE12 = zeros(4,modes);

for i = 1:modes
    FxE1(:,i) = K_column_x*Uxx1(:,i);
    FxE2(:,i) = K_beam_x*Uxx2(:,i);
    FxE3(:,i) = K_column_x*Uxx3(:,i);
    FxE4(:,i) = K_column_x*Uxx4(:,i);
    FxE5(:,i) = K_beam_x*Uxx5(:,i);
    FxE6(:,i) = K_column_x*Uxx6(:,i);
    FxE7(:,i) = K_column_x*Uxx7(:,i);
    FxE8(:,i) = K_beam_x*Uxx8(:,i);
    FxE9(:,i) = K_column_x*Uxx9(:,i);
    FxE10(:,i) = K_column_x*Uxx10(:,i);
    FxE11(:,i) = K_beam_x*Uxx11(:,i);
    FxE12(:,i) = K_column_x*Uxx12(:,i);
end
Forces_x=[FxE1;FxE2;FxE3;FxE4;FxE5;FxE6;FxE7;FxE8;FxE9;FxE10;FxE11;FxE12];
  
% y-axis
% The dispclacement related to the design spectrum:
Uyy1 = zeros(4,modes);
Uyy2 = zeros(4,modes);
Uyy3 = zeros(4,modes);
Uyy4 = zeros(4,modes);
Uyy5 = zeros(4,modes);
Uyy6 = zeros(4,modes);
Uyy7 = zeros(4,modes);
Uyy8 = zeros(4,modes);
Uyy9 = zeros(4,modes);
Uyy10 = zeros(4,modes);
Uyy11 = zeros(4,modes);
Uyy12 = zeros(4,modes);

for i = 1:modes
    Uyy1(:,i) = [u_d_y(1,i), u_d_y(3,i), 0, 0];
    Uyy2(:,i) = [0, u_d_y(3,i), 0, u_d_y(2,i)];
    Uyy3(:,i) = [u_d_y(1,i), u_d_y(2,i), 0, 0];
    Uyy4(:,i) = [u_d_y(4,i), u_d_y(6,i), u_d_y(1,i), u_d_y(3,i)];
    Uyy5(:,i) = [0, u_d_y(6,i), 0, u_d_y(5,i)];
    Uyy6(:,i) = [u_d_y(4,i), u_d_y(5,i), u_d_y(1,i), u_d_y(2,i)];
    Uyy7(:,i) = [u_d_y(7,i), u_d_y(9,i), u_d_y(4,i), u_d_y(6,i)];
    Uyy8(:,i) = [0, u_d_y(9,i), 0, u_d_y(8,i)];
    Uyy9(:,i) = [u_d_y(7,i), u_d_y(8,i), u_d_y(4,i), u_d_y(5,i)];
    Uyy10(:,i) = [u_d_y(10,i), u_d_y(12,i), u_d_y(7,i), u_d_y(9,i)];
    Uyy11(:,i) = [0, u_d_y(12,i), 0, u_d_y(11,i)];
    Uyy12(:,i) = [u_d_y(10,i), u_d_y(11,i), u_d_y(7,i), u_d_y(8,i)];
end


% ned fra og op
% for i = 1:modes
%     Uyy1(:,i) = [0, 0, u_d_y(1,i), u_d_y(3,i)];
%     Uyy2(:,i) = [0, u_d_y(3,i), 0, u_d_y(2,i)];
%     Uyy3(:,i) = [0, 0, u_d_y(1,i), u_d_y(2,i)];
%     Uyy4(:,i) = [u_d_y(1,i), u_d_y(3,i), u_d_y(4,i), u_d_y(6,i)];
%     Uyy5(:,i) = [0, u_d_y(6,i), 0, u_d_y(5,i)];
%     Uyy6(:,i) = [u_d_y(1,i), u_d_y(2,i), u_d_y(4,i), u_d_y(5,i)];
%     Uyy7(:,i) = [u_d_y(4,i), u_d_y(6,i), u_d_y(7,i), u_d_y(9,i)];
%     Uyy8(:,i) = [0, u_d_y(9,i), 0, u_d_y(8,i)];
%     Uyy9(:,i) = [u_d_y(4,i), u_d_y(5,i), u_d_y(7,i), u_d_y(8,i)];
%     Uyy10(:,i) = [u_d_y(7,i), u_d_y(9,i), u_d_y(10,i), u_d_y(12,i)];
%     Uyy11(:,i) = [0, u_d_y(12,i), 0, u_d_y(11,i)];
%     Uyy12(:,i) = [u_d_y(7,i), u_d_y(8,i), u_d_y(10,i), u_d_y(11,i)];
% end

% The local stiffness matriy for coloumn
K_column_y = E*Iy*[12/H^3,6/H^2,-12/H^3,6/H^2;
                  6/H^2,4/H,-6/H^2,2/H;
                  -12/H^3,-6/H^2,12/H^3,-6/H^2;
                  6/H^2,2/H,-6/H^2,4/H];

% The local stiffness matriy for beam
K_beam_y = E*Ib*[12/L^3,6/L^2,-12/L^3,6/L^2;
                  6/L^2,4/L,-6/L^2,2/L;
                  -12/L^3,-6/L^2,12/L^3,-6/L^2;
                  6/L^2,2/L,-6/L^2,4/L];


% Sectional forces for each element is found:
FyE1 = zeros(4,modes);
FyE2 = zeros(4,modes);
FyE3 = zeros(4,modes);
FyE4 = zeros(4,modes);
FyE5 = zeros(4,modes);
FyE6 = zeros(4,modes);
FyE7 = zeros(4,modes);
FyE8 = zeros(4,modes);
FyE9 = zeros(4,modes);
FyE10 = zeros(4,modes);
FyE11 = zeros(4,modes);
FyE12 = zeros(4,modes);

for i = 1:modes
    FyE1(:,i) = K_column_y*Uyy1(:,i);
    FyE2(:,i) = K_beam_y*Uyy2(:,i);
    FyE3(:,i) = K_column_y*Uyy3(:,i);
    FyE4(:,i) = K_column_y*Uyy4(:,i);
    FyE5(:,i) = K_beam_y*Uyy5(:,i);
    FyE6(:,i) = K_column_y*Uyy6(:,i);
    FyE7(:,i) = K_column_y*Uyy7(:,i);
    FyE8(:,i) = K_beam_y*Uyy8(:,i);
    FyE9(:,i) = K_column_y*Uyy9(:,i);
    FyE10(:,i) = K_column_y*Uyy10(:,i);
    FyE11(:,i) = K_beam_y*Uyy11(:,i);
    FyE12(:,i) = K_column_y*Uyy12(:,i);
end
Forces_y=[FyE1;FyE2;FyE3;FyE4;FyE5;FyE6;FyE7;FyE8;FyE9;FyE10;FyE11;FyE12];

%% FInd the forces with SRSS
modes = 12;
F_forcesSRSS_x1 = zeros(length(Forces_x),1);
F_forcesSRSS_x2 = zeros(length(Forces_x),1);
% Moments
for i = 1:length(Forces_x)
    for j = 1:modes
        F_forcesSRSS_x2(j) = Forces_x(i,j)^2;
        F_forcesSRSS_x1(i) = F_forcesSRSS_x1(i)+F_forcesSRSS_x2(j);
    end
end

F_forcesSRSS_x = sqrt(F_forcesSRSS_x1);

% y-axis
modes = 12;
F_forcesSRSS_y1 = zeros(length(Forces_y),1);
F_forcesSRSS_y2 = zeros(length(Forces_y),1);
% Moments
for i = 1:length(Forces_y)
    for j = 1:modes
        F_forcesSRSS_y2(j) = Forces_y(i,j)^2;
        F_forcesSRSS_y1(i) = F_forcesSRSS_y1(i)+F_forcesSRSS_y2(j);
    end
end

F_forcesSRSS_y = sqrt(F_forcesSRSS_y1);

%% FInd the forces with CQC
Beta_x = zeros(DOF,modes);
Beta_y = zeros(DOF,modes);
rho_x = zeros(DOF,modes);
rho_y = zeros(DOF,modes);
for i = 1:DOF
    for j=1:modes
        Beta_x(i,j)=omega_x(i)/omega_x(j);
        Beta_y(i,j)=omega_y(i)/omega_y(j);
        rho_x(i,j)=(8*(xi^2)*(1+Beta_x(i,j)*(Beta_x(i,j)^(3/2))))/((1-Beta_x(i,j)^2)^2+4*xi^2*Beta_x(i,j)*(1+Beta_x(i,j))^2);
        rho_y(i,j)=(8*(xi^2)*(1+Beta_y(i,j)*(Beta_y(i,j)^(3/2))))/((1-Beta_y(i,j)^2)^2+4*xi^2*Beta_y(i,j)*(1+Beta_y(i,j))^2);
    end
end
% x-axis 
modes = 12;
F_forceCQ1_x = zeros(length(Forces_x),1);
for i = 1:length(Forces_x)
    for j = 1:modes
        for n = 1:modes
            if n ~= j
                F_forceCQ1_x(i) = F_forceCQ1_x(i)+rho_x(j,n)*Forces_x(i,j)*Forces_x(i,n);
            else 
                F_forceCQ1_x(i) = F_forceCQ1_x(i);
            end
        end
    end
end


F_forceCQ_x = sqrt(F_forcesSRSS_x1+F_forceCQ1_x);

% For y-axis
modes = 12;
F_forceCQ1_y = zeros(length(Forces_y),1);
for i = 1:length(Forces_y)
    for j = 1:modes
        for n = 1:modes
            if n ~= j
                F_forceCQ1_y(i) = F_forceCQ1_y(i)+rho_y(j,n)*Forces_y(i,j)*Forces_y(i,n);
            else 
                F_forceCQ1_y(i) = F_forceCQ1_y(i);
            end
        end
    end
end


F_forceCQ_y = sqrt(F_forcesSRSS_y1+F_forceCQ1_y);

F_force_total = [F_forcesSRSS_x,F_forceCQ_x,F_forcesSRSS_y,F_forceCQ_y]


%% Exercise 3
TH2 = readmatrix('NORTHR_MU2125.txt');
TH = zeros(length(length(TH2)*5));
g = 9.81;
for j = 1:length(TH2)
      TH(1+(j-1)*5) = TH2(j,1);
      TH(2+(j-1)*5) = TH2(j,2);
      TH(3+(j-1)*5) = TH2(j,3);
      TH(4+(j-1)*5) = TH2(j,4);
      TH(5+(j-1)*5) = TH2(j,5);
end
TH = TH*g
x_TH = linspace(0,0.01*length(TH),length(TH));

PGA = max(abs(TH));

fig=figure()
fig.Position=[100 100 700 400];
plot(x_TH,TH(1,:))
hold on
plot(find(TH==-PGA)*0.01,-PGA,'r*')
xlabel('Time [s]')
hold off

% Scaling the plot
scaled = (TH .* ag *S)./PGA;

fig=figure()
fig.Position=[100 100 700 400];
hold on
plot(x_TH,TH(1,:))
plot(x_TH,scaled(1,:))
legend('Hector earthquake','Hector earthquake (scaled)')
hold off
%% Duhamel integral
modes = 12;
DOF = 12;
% Modal mass diag
M_n_x1 = diag(M_n_x);
M_n_y1 = diag(M_n_y);

% FInd delta tau which is the 1/20 of the lowest vibration period
Delta_tau_x = 1/20*min(Tx);
Delta_tau_y = 1/20*min(Ty);

% the damped natural angular frequency
omega_d_x = omega_x*sqrt(1-xi^2);
omega_d_y = omega_y*sqrt(1-xi^2);

% First thing to do is to find the constants M_1, M_2 and F
M_1x = 4*exp((-xi)*omega_x*Delta_tau_x);
M_2x = exp((-2)*xi*omega_x*Delta_tau_x);
Fx = Delta_tau_x./(3*M_n_x1.*omega_d_x);

M_1y = 4*exp((-xi)*omega_y*Delta_tau_y);
M_2y = exp((-2)*xi*omega_y*Delta_tau_y);
Fy = Delta_tau_y./(3*M_n_y1.*omega_d_y);

% Finding y_N and z_N for all the time step 
% First a time vector 
time_vector_x = 0:Delta_tau_x:x_TH(end);
time_vector_y = 0:Delta_tau_y:x_TH(end);

% Scale the data up
scaled_tme_vector_x = spline(x_TH,scaled,time_vector_x);
scaled_tme_vector_y = spline(x_TH,scaled,time_vector_y);


f_N_x = zeros(length(time_vector_x),length(M_n_x1));
f_N_y = zeros(length(time_vector_y),length(M_n_y1));
% Find the force f_N
for i = 1:length(time_vector_x)
    for j = 1:modes
        f_N_x(i,j) = scaled_tme_vector_x(i)*M_n_x1(j);
    end
end

for i = 1:length(time_vector_y)
    for j = 1:modes
        f_N_y(i,j) = scaled_tme_vector_y(i)*M_n_y1(j);
    end
end

% Finding y_N and z_N
y_N_x = zeros(length(time_vector_x),length(M_n_x1));
y_N_y = zeros(length(time_vector_y),length(M_n_y1));
z_N_x = zeros(length(time_vector_x),length(M_n_x1));
z_N_y = zeros(length(time_vector_y),length(M_n_y1));

for i = 1:length(scaled_tme_vector_x)
    for j = 1:modes
        y_N_x(i,j) = f_N_x(i,j) * cos(omega_d_x(j)*time_vector_x(i));
        z_N_x(i,j) = f_N_x(i,j) * sin(omega_d_x(j)*time_vector_x(i));
    end
end


for i = 1:length(scaled_tme_vector_y)
    for j = 1:modes
        y_N_y(i,j) = f_N_y(i,j) * cos(omega_d_y(j)*time_vector_y(i));
        z_N_y(i,j) = f_N_y(i,j) * sin(omega_d_y(j)*time_vector_y(i));
    end
end

% Now calculate A_N=2 and B_N=2
A_2_x = zeros(1,modes);
A_2_y = zeros(1,modes);
B_2_x = zeros(1,modes);
B_2_y = zeros(1,modes);

for i = 1:modes
    A_2_x(1,i) = Fx(i)*(y_N_x(1,i)*M_1x(i)+y_N_x(2,i));
    A_2_y(1,i) = Fy(i)*(y_N_y(1,i)*M_1y(i)+y_N_y(2,i));

    B_2_x(1,i) = Fx(i)*(z_N_x(1,i)*M_1x(i)+z_N_x(2,i));
    B_2_y(1,i) = Fy(i)*(z_N_y(1,i)*M_1y(i)+z_N_y(2,i));
end

% Now for A_N and B_N
A_N_x = zeros(length(scaled_tme_vector_x),modes);
A_N_y = zeros(length(scaled_tme_vector_y),modes);
B_N_x = zeros(length(scaled_tme_vector_x),modes);
B_N_y = zeros(length(scaled_tme_vector_y),modes);

A_N_x(2,:) = A_2_x;
A_N_y(2,:) = A_2_y;
B_N_x(2,:) = A_2_x;
B_N_y(2,:) = B_2_y;

for i = 4:2:length(scaled_tme_vector_x)
    for j = 1:modes
        A_N_x(i,j) = A_N_x(i-2,j)*M_2x(j)+Fx(j)*(y_N_x(i-2,j)*M_2x(j)+y_N_x(i-1,j)*M_2x(j)+y_N_x(i,j));
        B_N_x(i,j) = B_N_x(i-2,j)*M_2x(j)+Fx(j)*(z_N_x(i-2,j)*M_2x(j)+z_N_x(i-1,j)*M_1x(j)+z_N_x(i,j));
    end
end

for i = 4:2:length(scaled_tme_vector_y)
    for j = 1:modes
        A_N_y(i,j) = A_N_y(i-2,j)*M_2y(j)+Fy(j)*(y_N_y(i-2,j)*M_2y(j)+y_N_y(i-1,j)*M_1y(j)+y_N_y(i,j));
        B_N_y(i,j) = B_N_y(i-2,j)*M_2y(j)+Fy(j)*(z_N_y(i-2,j)*M_2y(j)+z_N_y(i-1,j)*M_1y(j)+z_N_y(i,j));
    end
end

% Now caclualte the response of a damped system
x_N_x = zeros(length(scaled_tme_vector_x),modes);
x_N_y = zeros(length(scaled_tme_vector_y),modes);

for i = 4:2:length(scaled_tme_vector_x)
    for j = 1:modes
        x_N_x(i,j) = A_N_x(i,j) * sin(omega_d_x(j)*time_vector_x(i))- B_N_x(i,j) * cos(omega_d_x(j)*time_vector_x(i));
    end
end

for i = 4:2:length(scaled_tme_vector_y)
    for j = 1:modes
        x_N_y(i,j) = A_N_y(i,j) * sin(omega_d_y(j)*time_vector_y(i)) - B_N_y(i,j) * cos(omega_d_y(j)*time_vector_y(i));
    end
end

plot(time_vector_y,x_N_y(:,1))

% All zero is removed from response of a damped system

x_N_x1 = zeros(length(scaled_tme_vector_x)/2,modes);
x_N_y1 = zeros(length(scaled_tme_vector_y)/2,modes);

for i = 1:length(scaled_tme_vector_x)/2-1
    for j = 1:modes
            x_N_x1(i,j) = x_N_x(i*2,j);
    end
end

for i = 1:length(scaled_tme_vector_y)/2-1
    for j = 1:modes
            x_N_y1(i,j) = x_N_y(i*2,j);
    end
end

time_vector_x1 = zeros(length(time_vector_x)/2,1);
time_vector_y1 = zeros(length(time_vector_y)/2,1);
% New time
for i = 1:(length(time_vector_x)/2)
      time_vector_x1(i) = time_vector_x(i*2);
end

for i = 1:(length(time_vector_y)/2)
      time_vector_y1(i) = time_vector_y(i*2);
end

figure ()
plot(time_vector_y1,x_N_y1(:,1))



% THe modeal displacements can be calculated u_n (n for each modes)
u_n_xDOF1 = zeros(length(x_N_x1),DOF);
u_n_xDOF2 = zeros(length(x_N_x1),DOF);
u_n_xDOF3 = zeros(length(x_N_x1),DOF);
u_n_xDOF4 = zeros(length(x_N_x1),DOF);
u_n_xDOF5 = zeros(length(x_N_x1),DOF);
u_n_xDOF6 = zeros(length(x_N_x1),DOF);
u_n_xDOF7 = zeros(length(x_N_x1),DOF);
u_n_xDOF8 = zeros(length(x_N_x1),DOF);
u_n_xDOF9 = zeros(length(x_N_x1),DOF);
u_n_xDOF10 = zeros(length(x_N_x1),DOF);
u_n_xDOF11 = zeros(length(x_N_x1),DOF);
u_n_xDOF12 = zeros(length(x_N_x1),DOF);
u_n_yDOF1= zeros(length(x_N_y1),DOF);
u_n_yDOF2 = zeros(length(x_N_y1),DOF);
u_n_yDOF3 = zeros(length(x_N_y1),DOF);
u_n_yDOF4= zeros(length(x_N_y1),DOF);
u_n_yDOF5 = zeros(length(x_N_y1),DOF);
u_n_yDOF6 = zeros(length(x_N_y1),DOF);
u_n_yDOF7= zeros(length(x_N_y1),DOF);
u_n_yDOF8 = zeros(length(x_N_y1),DOF);
u_n_yDOF9 = zeros(length(x_N_y1),DOF);
u_n_yDOF10= zeros(length(x_N_y1),DOF);
u_n_yDOF11 = zeros(length(x_N_y1),DOF);
u_n_yDOF12 = zeros(length(x_N_y1),DOF);

for i = 1:modes
    u_n_xDOF1(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(1,i));
    u_n_xDOF2(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(2,i));
    u_n_xDOF3(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(3,i));
    u_n_xDOF4(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(4,i));
    u_n_xDOF5(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(5,i));
    u_n_xDOF6(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(6,i));
    u_n_xDOF7(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(7,i));
    u_n_xDOF8(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(8,i));
    u_n_xDOF9(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(9,i));
    u_n_xDOF10(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(10,i));
    u_n_xDOF11(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(11,i));
    u_n_xDOF12(:,i) = (x_N_x1(:,i)*Gamma_n_x(i)*phi_x(12,i));
end

for i = 1:modes
    u_n_yDOF1(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(1,i);
    u_n_yDOF2(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(2,i);
    u_n_yDOF3(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(3,i);
    u_n_yDOF4(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(4,i);
    u_n_yDOF5(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(5,i);
    u_n_yDOF6(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(6,i);
    u_n_yDOF7(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(7,i);
    u_n_yDOF8(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(8,i);
    u_n_yDOF9(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(9,i);
    u_n_yDOF10(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(10,i);
    u_n_yDOF11(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(11,i);
    u_n_yDOF12(:,i) = x_N_y1(:,i)*Gamma_n_y(i)*phi_y(12,i);
end



% Then the total displacement is found U for DOF for each DOF
U_xDOF1 = nan(length(x_N_x1),1);
U_xDOF2 = nan(length(x_N_x1),1);
U_xDOF3 = nan(length(x_N_x1),1);
U_xDOF4 = nan(length(x_N_x1),1);
U_xDOF5 = nan(length(x_N_x1),1);
U_xDOF6 = nan(length(x_N_x1),1);
U_xDOF7 = nan(length(x_N_x1),1);
U_xDOF8 = nan(length(x_N_x1),1);
U_xDOF9 = nan(length(x_N_x1),1);
U_xDOF10 = nan(length(x_N_x1),1);
U_xDOF11 = nan(length(x_N_x1),1);
U_xDOF12 = nan(length(x_N_x1),1);


U_yDOF1 = nan(length(x_N_y1),1);
U_yDOF2 = nan(length(x_N_y1),1);
U_yDOF3 = nan(length(x_N_y1),1);
U_yDOF4 = nan(length(x_N_y1),1);
U_yDOF5 = nan(length(x_N_y1),1);
U_yDOF6 = nan(length(x_N_y1),1);
U_yDOF7 = nan(length(x_N_y1),1);
U_yDOF8 = nan(length(x_N_y1),1);
U_yDOF9 = nan(length(x_N_y1),1);
U_yDOF10 = nan(length(x_N_y1),1);
U_yDOF11 = nan(length(x_N_y1),1);
U_yDOF12 = nan(length(x_N_y1),1);

for i = 1:length(x_N_x1)
    U_xDOF1(i) = sum(u_n_xDOF1(i,:));
    U_xDOF2(i) = sum(u_n_xDOF2(i,:));
    U_xDOF3(i) = sum(u_n_xDOF3(i,:));
    U_xDOF4(i) = sum(u_n_xDOF4(i,:));
    U_xDOF5(i) = sum(u_n_xDOF5(i,:));
    U_xDOF6(i) = sum(u_n_xDOF6(i,:));
    U_xDOF7(i) = sum(u_n_xDOF7(i,:));
    U_xDOF8(i) = sum(u_n_xDOF8(i,:));
    U_xDOF9(i) = sum(u_n_xDOF9(i,:));
    U_xDOF10(i) = sum(u_n_xDOF10(i,:));
    U_xDOF11(i) = sum(u_n_xDOF11(i,:));
    U_xDOF12(i) = sum(u_n_xDOF12(i,:));
end

U_x = [U_xDOF1,U_xDOF2,U_xDOF3,U_xDOF4,U_xDOF5,U_xDOF6,U_xDOF7,U_xDOF8,U_xDOF9,U_xDOF10,U_xDOF11,U_xDOF12];

for i = 1:length(x_N_y1)
    U_yDOF1(i) = sum(u_n_yDOF1(i,:));
    U_yDOF2(i) = sum(u_n_yDOF2(i,:));
    U_yDOF3(i) = sum(u_n_yDOF3(i,:));
    U_yDOF4(i) = sum(u_n_yDOF4(i,:));
    U_yDOF5(i) = sum(u_n_yDOF5(i,:));
    U_yDOF6(i) = sum(u_n_yDOF6(i,:));
    U_yDOF7(i) = sum(u_n_yDOF7(i,:));
    U_yDOF8(i) = sum(u_n_yDOF8(i,:));
    U_yDOF9(i) = sum(u_n_yDOF9(i,:));
    U_yDOF10(i) = sum(u_n_yDOF10(i,:));
    U_yDOF11(i) = sum(u_n_yDOF11(i,:));
    U_yDOF12(i) = sum(u_n_yDOF12(i,:));
end
U_y = [U_yDOF1,U_yDOF2,U_yDOF3,U_yDOF4,U_yDOF5,U_yDOF6,U_yDOF7,U_yDOF8,U_yDOF9,U_yDOF10,U_yDOF11,U_yDOF12];


%%% PLOTTING
% x direction
fig = figure ();
fig.Position=[100 100 1600 700]
for i = 1:DOF
    subplot(3,4,i)
    hold on
    plot(time_vector_x1,U_x(:,i))
    plot(time_vector_x1(find(U_x(:,i)==max(U_x(:,i)),1)),max(U_x(:,i)),'xr')
    if (i==1 ||  i==4 ||  i==7 || i== 10)
        ylabel('Displacement [m]')
    else
        ylabel('Rotation [\circ]')
    end
    xlabel('Time [s]')
    yline(r_tot_xe(i),'--r')
    yline(-r_tot_xe(i),'--r')
    title(sprintf('DOF %d',i))
    %axis padded
end
% han=axes(fig,'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'Displacement [m]');
% xlabel(han,'Time [s]');
saveas(gcf,'C:\Users\Frede\Danmarks Tekniske Universitet\Katarina Kolding Skaarup - 41967 Seismic and wind engineering\Figures\Duhamel_x.eps','epsc')

% y direction
fig1 = figure ();
fig1.Position=[100 100 1600 700]
for i = 1:DOF
    subplot(3,4,i)
    hold on
    plot(time_vector_y1,U_y(:,i))
    plot(time_vector_y1(find(U_y(:,i)==max(U_y(:,i)),1)),max(U_y(:,i)),'xr')
    yline(r_tot_ye(i),'--r')
    yline(-r_tot_ye(i),'--r')
    if (i==1 || i==4 || i==7 || i== 10)
        ylabel('Displacement [m]')
    else
        ylabel('Rotation [\circ]')
    end
    xlabel('Time [s]')
    title(sprintf('DOF %d',i))
    %axis padded
end
% han=axes(fig1,'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'Displacement [m]');
% xlabel(han,'Time [s]');
saveas(gcf,'C:\Users\Frede\Danmarks Tekniske Universitet\Katarina Kolding Skaarup - 41967 Seismic and wind engineering\Figures\Duhamel_y.eps','epsc')


% difference in deformation
et = zeros(12,1);
to = zeros(12,1);
tre = zeros(12,1);
fire = zeros(12,1);
diffx = zeros(12,1);
diffy = zeros(12,1);
for i=1:12
et(i)=max(r_tot_xeCQ(i));
to(i) = max(U_x(:,i));
diffx(i) = (et(i)-to(i))/to(i)*100;
end

for i=1:12
tre(i)=max(r_tot_yeCQ(i));
fire(i) = max(U_y(:,i));
diffy(i) = (tre(i)-fire(i))/fire(i)*100;
end
%% Now we want to compare MRSA and THA
scaled_tme_vector_x;   % Load montion, First Col. is the time and the second col. is Acce. in g
g=9.810;             % (m/s^2)
Ag=scaled*g; % Acceleration Vector
dt=0.01;            % Time Interval (sec)
zet=xi*100;              % Damping Ratio (%)
endp=4;             % End Period of Spectra (sec) 

[T_new,Spa_new,Spv_new,Sd_new]=SPEC(dt,Ag,zet,g,endp);

% figure ()
% plot(T_new,Spa_new)
% grid on
% xlabel('Period (sec)','FontSize',13);
% ylabel('Spa (g)','FontSize',13);
% title('Displacement Spectrum','FontSize',13)

% fig = figure()
% fig.Position=[100 100 1600 700]
figure ()
plot(T_new,Spa_new,'-r')
hold on
plot(T,Se,'-b')
xline(Tx(1),'--') % x-direction
xline(Ty(1),'-.') % y-direction
hold off
grid on
legend('THA','MSRS','First period x-direction','First period y-direction')
xlabel('Period, T [s]')
ylabel('Spectral Acceleration S_e [g]')
saveas(gcf,'C:\Users\Frede\Danmarks Tekniske Universitet\Katarina Kolding Skaarup - 41967 Seismic and wind engineering\Figures\compareMRSATSA_x.eps','epsc')

%% Exercise 4 - Modal drift 
drift_x = zeros(4,12);
drift_y = zeros(4,12);

for i = 1:12 
    drift_x(:,i) = [u_d_x(1,i); u_d_x(4,i)-u_d_x(1,i); u_d_x(7,i)-u_d_x(4,i); u_d_x(10,i)-u_d_x(7,i)];
    drift_y(:,i) = [u_d_y(1,i); u_d_y(4,i)-u_d_y(1,i); u_d_y(7,i)-u_d_y(4,i); u_d_y(10,i)-u_d_y(7,i)] ;
end

% x-axis 
DOF = 12;
modes = 12;

drift_SRSS_x1 = zeros(4,1);
drift_SRSS_x2 = zeros(modes,1);
drift_SRSS_y1 = zeros(4,1);
drift_SRSS_y2 = zeros(modes,1);
for i = 1:4
    for j = 1:modes
        drift_SRSS_x2(j) = drift_x(i,j)^2;
        drift_SRSS_x1(i) = drift_SRSS_x1(i)+drift_SRSS_x2(j);
        drift_SRSS_y2(j) = drift_y(i,j)^2;
        drift_SRSS_y1(i) = drift_SRSS_y1(i)+drift_SRSS_y2(j);
    end
end
drift_SRSS_x = sqrt(drift_SRSS_x1);
drift_SRSS_y = sqrt(drift_SRSS_y1);

drift_CQC_x1 = zeros(4,1);
drift_CQC_y1 = zeros(4,1);
for i = 1:4
    for j = 1:modes
        for n = 1:modes
            if n ~= j
                drift_CQC_x1(i) = drift_CQC_x1(i)+rho_x(j,n)*drift_x(i,j)*drift_x(i,n);
                drift_CQC_y1(i) = drift_CQC_y1(i)+rho_y(j,n)*drift_y(i,j)*drift_y(i,n);
            else 
                drift_CQC_x1(i) = drift_CQC_x1(i);
                drift_CQC_y1(i) = drift_CQC_y1(i);
            end
        end
    end
end

drift_QCQ_x = sqrt(drift_SRSS_x1+drift_CQC_x1)*q*1000; % [mm]
drift_QCQ_y = sqrt(drift_SRSS_y1+drift_CQC_y1)*q*1000; % [mm]

max_drift = 0.0075*H/0.5; %i m

%u_d_yend