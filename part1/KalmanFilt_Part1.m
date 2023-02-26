clear; % Clear variables
datasetNum = 9; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime] = init(datasetNum);
Z = sampledVicon(1:6,:);%all the measurements that you need for the update
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %J ust for saving state his.
prevTime = 0; %last time step in real time
%write your code here calling the pred_step.m and upd_step.m functions
%Name: Pavan Chowdary Cherukuri
%NetId: pc3088
%N number: N10938396
%Kalman Filter Part 1
for i = 1:length(sampledTime)


    %Prediction step

    %Inputs-

    %uPrev is given by uPrev

    %covarPrev is given by covarPrev

    %angVel is given by sampledData(i).omg which represents the ith column
    %is extraced from omg values in sampledData

    %acc is given by sampledData(i).acc which represents the ith column is
    %extracted from acc values in sampledData

    %dt is given by sampledTime(i)-prevTime

    [covarEst,uEst] = pred_step(uPrev,covarPrev,sampledData(i).omg,sampledData(i).acc,sampledTime(i)-prevTime);
    


    %Update Step

    %Inputs-

    %z_t is given by sampledVicon(1:6,i) which is the ith column of
    %velocity components from Vicon

    %covarEst is given by covarEst

    %uEst is given by uEst

    [uCurr,covar_curr] = upd_step(sampledVicon(1:6,i),covarEst,uEst);
    



    %prevTime is re initialised to sampleTime(i) for calculation of dt in 
    %next iteration

    prevTime = sampledTime(i);

    %The calculated uCurr from the update step is saved to savedStates ith
    %column

    savedStates(:,i) = uCurr;

    %uPrev and covarPrev are reinitialised as uCurr and covar_curr 
    %for next iteration

    uPrev = uCurr;
    covarPrev = covar_curr;


end



plotData(savedStates, sampledTime, sampledVicon, 1, datasetNum);