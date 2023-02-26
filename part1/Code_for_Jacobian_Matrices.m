%Name: Pavan Chowdary Cherukuri
%NetId: pc3088
%N number: N10938396

%Code for calculating the Jacobian Matrices

%symbolic variables roll,pitch,yaw were defined
syms roll
syms pitch
syms yaw

%Rotation Matrices were defined
RY = roty(pitch) 
RZ = rotz(yaw) 
RX = rotx(roll)

%Calculation of the RZYX- R(x2) matrix
RZYX = RZ*RY*RX; 

%gravity vector
g = [0; 0; -9.8]; 

%G(x2) Matrix
GZYX = [ cos(pitch)*cos(yaw) -sin(yaw) 0;
         sin(yaw)*cos(pitch) cos(yaw) 0;
         -sin(pitch)  0  1]

%inverse of G(x2) Matrix
GZYX_inv = inv(GZYX); 

%Now, the process model 'x' is given as-
% [x1;x2;x3;x4;x5] matrix, where,x1 = p,x2 = q,x3 = p_dot,x4 = bg,x5 = ba
%'p' is the position, 'q' is the orientation, 'p_dot' is the linear
%velocity, 'bg' is the gyroscope bias, 'ba' is the accelerometer bias.
% 'x' is a 15*1 matrix

%The Process model will be given by x_dot = [x3;
                                            %GZYX_inv*RZYX*(wm-x4-ng);
                                            %g+R(x2)*(am-x5-na);
                                            %nbg;
                                            %nba] = f(x,u,n)
%nbg and nba are the drift in the gyroscope and accelerometer bias 
%x_dot is given by function 'f' which depends on x,u,n; where 'x' is the
%state vector, 'u' is control inputs- angular velocity and acceleration
%'n' is the matrix that consists of bias values- [0;
                                                 %ng;
                                                 %na;
                                                 %nbg;
                                                 %nba]


%Defining various states and bias values as symbolic variables for
%calculating the Jacobian Matrix

x1 = sym('x1',[3 1]); 
x2 = [roll; pitch; yaw];
vm = sym('vm',[3 1]); 
bg = sym('bg',[3 1]);
ba = sym('ba',[3 1]);
wm = sym('wm',[3 1]); 
am = sym('am',[3 1]);
ng = sym('ng',[3 1]);
na = sym('na',[3 1]);
nbg = sym('nbg',[3 1]);
nba = sym('nba',[3 1]);

%%
%Calculating dx_dot/dx or df/dx - Jacobian matrix:

%Here the representation of the non zero elements in the Jacobian Matrix
%were given by At_ij (i- row number, j- column number)
%It can be understood that all the zero elements indicate that the partial
%differentiation in the of 'f' w.r.t 'x' was zero for that element

%For finding At_44,At_54,At_64 below shown operation is performed which 
%is differentiating GZYX_inv*RZYX with roll and multiplying with (wm-bg)
At_i4 = (simplify(diff(simplify(GZYX_inv*RZYX),roll)))*(wm-bg)

At_44 = At_i4(1);
At_54 = At_i4(2);
At_64 = At_i4(3);

%For finding At_45,At_55,At_65 the below shown operation is performed which 
%is differentiating GZYX_inv*RZYX with pitch and multiplying with (wm-bg)
At_i5 = (simplify(diff(simplify(GZYX_inv*RZYX),pitch)))*(wm-bg)
At_45 = At_i5(1);
At_55 = At_i5(2);
At_65 = At_i5(3);

%For finding At_46,At_56,At_66 the below shown operation is performed which 
%is differentiating GZYX_inv*RZYX with yaw and multiplying with (wm-bg)
At_i6 = (simplify(diff(simplify(GZYX_inv*RZYX),yaw)))*(wm-bg)
At_46 = At_i6(1);
At_56 = At_i6(2);
At_66 = At_i6(3);

%For finding At_74,At_84,At_94 the below shown operation is performed which 
%is differentiating RZYX with roll and multiplying with (am-ba)
At_ia4 = (simplify(diff(RZYX,roll)))*(am-ba)
At_74 = At_ia4(1);
At_84 = At_ia4(2);
At_94 = At_ia4(3);

%For finding At_75,At_85,At_95 the below shown operation is performed which 
%is differentiating RZYX with pitch and multiplying with (am-ba)
At_ia5 = (simplify(diff(RZYX,pitch)))*(am-ba)
At_75 = At_ia5(1);
At_85 = At_ia5(2);
At_95 = At_ia5(3);

%For finding At_76,At_86,At_96 the below shown operation is performed which 
%is differentiating RZYX with yaw and multiplying with (am-ba)
At_ia6 = (simplify(diff(RZYX,yaw)))*(am-ba)
At_76 = At_ia6(1);
At_86 = At_ia6(2);
At_96 = At_ia6(3);

%When the first element in 'x_dot' which is 'x3' is differentiated with 
% respect to x3, then the value will be 1. Hence the elements At_17,
%At_28,At_39 will be equal to '1'


%For finding At_410,At_411,At_412,At_510,At_511,At_512,At_610,At_611,At_612
%differentiate inv(G(x2))*R(x2)*(wm-x4-ng) with respect to x4 
%which gives inv(G(x2))*R(x2)*(-I), which is shown below
At_x2dot_diff_x4 = GZYX_inv*RZYX*[-1,0,0;
                         0,-1,0;
                         0,0,-1]
At_410 = At_x2dot_diff_x4(1,1);
At_411 = At_x2dot_diff_x4(1,2);
At_412 = At_x2dot_diff_x4(1,3);
At_510 = At_x2dot_diff_x4(2,1);
At_511 = At_x2dot_diff_x4(2,2);
At_512 = At_x2dot_diff_x4(2,3);
At_610 = At_x2dot_diff_x4(3,1);
At_611 = At_x2dot_diff_x4(3,2);
At_612 = At_x2dot_diff_x4(3,3);

%For finding At_713,At_714,At_715,At_813,At_814,At_815,At_913,At_914,At_915
%differentiate R(x2)*(am-x5-na) with respect to x5 
%which gives R(x2)*(-I), which is shown below
At_x3dot_diff_x5 = RZYX*[-1,0,0;
                         0,-1,0;
                         0,0,-1]
At_713 = At_x3dot_diff_x5(1,1);
At_714 = At_x3dot_diff_x5(1,2);
At_715 = At_x3dot_diff_x5(1,3);
At_813 = At_x3dot_diff_x5(2,1);
At_814 = At_x3dot_diff_x5(2,2);
At_815 = At_x3dot_diff_x5(2,3);
At_913 = At_x3dot_diff_x5(3,1);
At_914 = At_x3dot_diff_x5(3,2);
At_915 = At_x3dot_diff_x5(3,3);


%%
%For calculating dfdn Jacobian Matrix:
%The [(4,4),(4,5),(4,6);
%     (5,4),(5,5),(5,6);
%     (6,4),(6,5),(6,6)] elements in dfdn are same as
% [At_410,At_411,At_412;
%  At_510,At_511,At_512;
%  At_610,At_611,At_612] elements in the At matrix
%and the [(7,7),(7,8),(7,9);
%         (8,7),(8,8),(8,9);
%         (9,7),(9,8),(9,9)] elements in dfdn are same as
% [At_713,At_714,At_715;
%  At_813,At_814,At_815;
%  At_913,At_914,At_915] elements in the At matrix
%This deduction can be easily understood as the partial derivatives
%these elements represent were the same.
%Hence these elements are directly given in the dfdn matrix
%The (10,10),(11,11),(12,12),(13,13),(14,14),(15,15) elements in dfdn are
%'1'

%%
%Part 1
%For Calculating Ct Jacobian matrix
%C matrix is given as
C =[ 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; 
     0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
    ];
%x matrix
x = [x1;
     x2;
     vm;
     bg;
     ba];
%z
z = C*x; %Since noise is constant when partially differentiated with 'x'
         %it is not considered here.
         %Actual Measurement model will be z = C*x+v
Ct = jacobian(z,x)




%%
%Part 2
%For Calculating Ct Jacobian matrix
%C matrix is given as
C =[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 
      0 0 0 0 0 0 0 0 1 0 0 0 0 0 0
    ];
%x matrix
x = [x1;
     x2;
     vm;
     bg;
     ba];
%z
z = C*x; %Since noise is constant when partially differentiated with 'x'
         %it is not considered here.
         %Actual Measurement model will be z = C*x+v
Ct = jacobian(z,x)









