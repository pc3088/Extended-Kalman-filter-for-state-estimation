 function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%covarPrev and uPrev are the previous mean and covariance respectively
%angVel is the angular velocity
%acc is the acceleration
%dt is the sampling time

wm1 = angVel(1,1);
wm2 = angVel(2,1);
wm3 = angVel(3,1);
wm = [wm1;wm2;wm3];
am1 = acc(1,1);
am2 = acc(2,1);
am3 = acc(3,1);
am = [am1;am2;am3];
x1 = uPrev(1,1);
y1 = uPrev(2,1);
z1 = uPrev(3,1);
roll= uPrev(4,1);
pitch= uPrev(5,1);
yaw= uPrev(6,1);
bg = uPrev(10:12);
bg1 = bg(1);
bg2 = bg(2);
bg3 = bg(3);
ba = uPrev(13:15);
ba1 = ba(1);
ba2 = ba(2);
ba3 = ba(3);
RY = roty(pitch); 
RZ = rotz(yaw); 
RX = rotx(roll);
RZYX = RZ*RY*RX;
g = [0; 0; -9.81];
GZYX = [ cos(pitch)*cos(yaw) -sin(yaw) 0;
         sin(yaw)*cos(pitch) cos(yaw)  0;
         -sin(pitch)         0         1];
GZYX_inv = inv(GZYX);

At_i4 = [
((sin(pitch)*sin(roll)*(bg3 - wm3))/cos(pitch) - (cos(roll)*sin(pitch)*(bg2 - wm2))/cos(pitch));
                                                (cos(roll)*(bg3 - wm3) + sin(roll)*(bg2 - wm2));
                      ((sin(roll)*(bg3 - wm3))/cos(pitch) - (cos(roll)*(bg2 - wm2))/cos(pitch))
 ];

At_44 = At_i4(1);
At_54 = At_i4(2);
At_64 = At_i4(3);

At_i5 =                [   ((sin(roll)*(bg2 - wm2))/(sin(pitch)^2 - 1) - (cos(roll)*(bg3 - wm3))/cos(pitch)^2);
                                                                                                      0;
((sin(pitch)*sin(roll)*(bg2 - wm2))/(sin(pitch)^2 - 1) - (cos(roll)*sin(pitch)*(bg3 - wm3))/cos(pitch)^2)
 ];
At_45 = At_i5(1);
At_55 = At_i5(2);
At_65 = At_i5(3);

At_i6 = [0;0;0];
At_46 = At_i6(1);
At_56 = At_i6(2);
At_66 = At_i6(3);

At_ia4 = [
  ((am2 - ba2)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)) + (am3 - ba3)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)));
(- (am2 - ba2)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)) - (am3 - ba3)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)));
                                                                  (cos(pitch)*cos(roll)*(am2 - ba2) - cos(pitch)*sin(roll)*(am3 - ba3))
                                                                  ];
At_74 = At_ia4(1);
At_84 = At_ia4(2);
At_94 = At_ia4(3);

At_ia5 = [
 
(cos(pitch)*cos(roll)*cos(yaw)*(am3 - ba3) - cos(yaw)*sin(pitch)*(am1 - ba1) + cos(pitch)*cos(yaw)*sin(roll)*(am2 - ba2));
(cos(pitch)*cos(roll)*sin(yaw)*(am3 - ba3) - sin(pitch)*sin(yaw)*(am1 - ba1) + cos(pitch)*sin(roll)*sin(yaw)*(am2 - ba2));
                         (- cos(pitch)*(am1 - ba1) - cos(roll)*sin(pitch)*(am3 - ba3) - sin(pitch)*sin(roll)*(am2 - ba2))
 ];
At_75 = At_ia5(1);
At_85 = At_ia5(2);
At_95 = At_ia5(3);

At_ia6 = [
 
((am3 - ba3)*(cos(yaw)*sin(roll) - cos(roll)*sin(pitch)*sin(yaw)) - (am2 - ba2)*(cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw)) - cos(pitch)*sin(yaw)*(am1 - ba1));
((am3 - ba3)*(sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch)) - (am2 - ba2)*(cos(roll)*sin(yaw) - cos(yaw)*sin(pitch)*sin(roll)) + cos(pitch)*cos(yaw)*(am1 - ba1));
                                                                                                                                                                    0
 ];
At_76 = At_ia6(1);
At_86 = At_ia6(2);
At_96 = At_ia6(3);

At_x2dot_diff_x4 =  [-cos(pitch)*cos(yaw),   cos(roll)*sin(yaw)-cos(yaw)*sin(pitch)*sin(roll), -sin(roll)*sin(yaw)-cos(roll)*cos(yaw)*sin(pitch);
-cos(pitch)*sin(yaw), -cos(roll)*cos(yaw)-sin(pitch)*sin(roll)*sin(yaw),   cos(yaw)*sin(roll)-cos(roll)*sin(pitch)*sin(yaw);
          sin(pitch),                                -cos(pitch)*sin(roll),                                -cos(pitch)*cos(roll)];
At_410 = At_x2dot_diff_x4(1,1);
At_411 = At_x2dot_diff_x4(1,2);
At_412 = At_x2dot_diff_x4(1,3);
At_510 = At_x2dot_diff_x4(2,1);
At_511 = At_x2dot_diff_x4(2,2);
At_512 = At_x2dot_diff_x4(2,3);
At_610 = At_x2dot_diff_x4(3,1);
At_611 = At_x2dot_diff_x4(3,2);
At_612 = At_x2dot_diff_x4(3,3);


At_x3dot_diff_x5 = [-cos(pitch)*cos(yaw),   cos(roll)*sin(yaw)-cos(yaw)*sin(pitch)*sin(roll), -sin(roll)*sin(yaw)-cos(roll)*cos(yaw)*sin(pitch);
-cos(pitch)*sin(yaw), -cos(roll)*cos(yaw)-sin(pitch)*sin(roll)*sin(yaw),   cos(yaw)*sin(roll)-cos(roll)*sin(pitch)*sin(yaw);
          sin(pitch),                                -cos(pitch)*sin(roll),                                -cos(pitch)*cos(roll)];
At_713 = At_x3dot_diff_x5(1,1);
At_714 = At_x3dot_diff_x5(1,2);
At_715 = At_x3dot_diff_x5(1,3);
At_813 = At_x3dot_diff_x5(2,1);
At_814 = At_x3dot_diff_x5(2,2);
At_815 = At_x3dot_diff_x5(2,3);
At_913 = At_x3dot_diff_x5(3,1);
At_914 = At_x3dot_diff_x5(3,2);
At_915 = At_x3dot_diff_x5(3,3);



dfdx =  [0,0,0,0,    0,    0,    1,0,0,0,     0,     0,     0,     0,     0;
         0,0,0,0,    0,    0,    0,1,0,0,     0,     0,     0,     0,     0;
         0,0,0,0,    0,    0,    0,0,1,0,     0,     0,     0,     0,     0;
         0,0,0,At_44,At_45,At_46,0,0,0,At_410,At_411,At_412,0,     0,     0;
         0,0,0,At_54,At_55,At_56,0,0,0,At_510,At_511,At_512,0,     0,     0;
         0,0,0,At_64,At_65,At_66,0,0,0,At_610,At_611,At_612,0,     0,     0;
         0,0,0,At_74,At_75,At_76,0,0,0,0,     0,     0,     At_713,At_714,At_715;
         0,0,0,At_84,At_85,At_86,0,0,0,0,     0,     0,     At_813,At_814,At_815;
         0,0,0,At_94,At_95,At_96,0,0,0,0,     0,     0,     At_913,At_914,At_915;
         0,0,0,0,    0,    0,    0,0,0,0,     0,     0,     0,     0,     0;
         0,0,0,0,    0,    0,    0,0,0,0,     0,     0,     0,     0,     0;
         0,0,0,0,    0,    0,    0,0,0,0,     0,     0,     0,     0,     0;
         0,0,0,0,    0,    0,    0,0,0,0,     0,     0,     0,     0,     0;
         0,0,0,0,    0,    0,    0,0,0,0,     0,     0,     0,     0,     0;
         0,0,0,0,    0,    0,    0,0,0,0,     0,     0,     0,     0,     0];
    


At = dfdx;


dfdn = [0,0,0,0,     0,     0,     0,     0,     0,     0,0,0,0,0,0;
        0,0,0,0,     0,     0,     0,     0,     0,     0,0,0,0,0,0;
        0,0,0,0,     0,     0,     0,     0,     0,     0,0,0,0,0,0;
        0,0,0,At_410,At_411,At_412,0,     0,     0,     0,0,0,0,0,0;
        0,0,0,At_510,At_511,At_512,0,     0,     0,     0,0,0,0,0,0;
        0,0,0,At_610,At_611,At_612,0,     0,     0,     0,0,0,0,0,0;
        0,0,0,0,     0,     0,     At_713,At_714,At_715,0,0,0,0,0,0;
        0,0,0,0,     0,     0,     At_813,At_814,At_815,0,0,0,0,0,0;
        0,0,0,0,     0,     0,     At_913,At_914,At_915,0,0,0,0,0,0;
        0,0,0,0,     0,     0,     0,     0,     0,     1,0,0,0,0,0;
        0,0,0,0,     0,     0,     0,     0,     0,     0,1,0,0,0,0;
        0,0,0,0,     0,     0,     0,     0,     0,     0,0,1,0,0,0;
        0,0,0,0,     0,     0,     0,     0,     0,     0,0,0,1,0,0;
        0,0,0,0,     0,     0,     0,     0,     0,     0,0,0,0,1,0;
        0,0,0,0,     0,     0,     0,     0,     0,     0,0,0,0,0,1];
        

Ut = dfdn;

Ft = eye(15)+ dt*At;
Vt = Ut;
Q = eye(15);
Qd = Q*dt;

vm = uPrev(7:9,1);
nbg = [0; 0; 0];
nba = [0; 0; 0];

fxun =  [vm;
         GZYX_inv*RZYX*(wm-bg-nbg); 
         g + RZYX*(am-ba-nba);
         nbg;
         nba];


uEst = uPrev + dt*fxun;
covarEst = Ft*covarPrev*(transpose(Ft)) + Vt*Qd*(transpose(Vt));




end

