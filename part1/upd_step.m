function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)


%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively
%Name: Pavan Chowdary Cherukuri
%NetId: pc3088
%N number: N10938396
%Update Step Part 1

%Process Model zt = Cx+v
%C = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%     0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; 
%     0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
%     0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 
%     0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
%     0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
%    ];
%since only the velocity is considered for the measurement model in
% update step in part 2

%Ct- Jacobian Matrix taken from "Code_for_Jacobian_Matrices.m" file


Ct = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; 
     0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 
     0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
    ];
%g_uEst is g(uEst,0) which in this case will be 7,8,9 elements of uEst
g_uEst = uEst(1:6,1);

%Matrix R is the co-variance matrix of noise which is given by 0.01*I
%Note- R is subject to change varying the values of the constant multiplied
%can vary different parameters
R = eye(6)*0.01;

%transpose of Ct is defined
Ct_transpose = transpose(Ct);

%Kt is defined
%In the formula 'Wt' and 'transpose(Wt)' are ignored as they are Identity
%matrices
Kt = covarEst*(Ct_transpose)*inv((Ct*covarEst*Ct_transpose)+ R);

%uCurr and covar_curr are calculated
uCurr = uEst + Kt*(z_t - g_uEst);
covar_curr = covarEst - Kt*Ct*covarEst;

end