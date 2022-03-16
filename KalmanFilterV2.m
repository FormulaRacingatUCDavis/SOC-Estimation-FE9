
clc;
%formatlatex
format longG
warning('off',  'all')
mkdir('./Figures')
warning('on',  'all')


%--------------------------------------------------------------------------
% At time k-1
%--------------------------------------------------------------------------

%
% Constants
% k = (time value);
% T = (value);
dt = 1;   %Sampling Period
R0 = 0.01;
Rc = 0.015;
Ccap = 2400;
Cbat = 18000;
Voc0 = 3.435;
alp = 0.65;

%
% Setting up test data
%Samples = 100;
%actualSOC = ones(1, Samples);
%V = ones(1, Samples);
%I = ones(1, Samples);
%timeSteps = ones(1, Samples);
% importing from csv
%T = csvread("data_sim.csv");
%for i = 1:Samples
%    timeSteps(i) = T(i,1);
%    actualSOC(i) = T(i,2);
%    V(i) = T(i,3);
%    I(i) = T(i,4);
%end

actualSOC = VarName2.';  %[0.3, 0.6, 0.5, 0.4, 0.5, 0.3, 0.6, 0.4]  %Actual SOC
Vc = VarName3.';
V = VarName4.';                          %[0.3, 0.6, 0.5, 0.4, 0.5, 0.3, 0.6, 0.4];         %Measured Voltage
I = VarName5.';                          %Measured Current
timeSteps = VarName1.';                  %Time
%totalTime = length(timeSteps);

%V(20) = 0;
%V(21) = 0;
%V(22) = 0;
%V(23) = 0;
%V(24) = 0;


%
%Simulation Data goes here
%
%actualSOC = [0.3, 0.6, 0.5, 0.4, 0.5, 0.3, 0.6, 0.4]  %Actual SOC
%V = [0.3, 0.6, 0.5, 0.4, 0.5, 0.3, 0.6, 0.4];         %Measured Voltage
%I = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];           %Measured Current
%timeSteps = [0:0.1:length(I)];                            %Time
totalTime = length(timeSteps);

%xhatk_1 = [SOC; Vc];
%I is measured

%
% Initializing probability matricies
%
Aprime = [1, 0; 0, exp(-dt/(Ccap*Rc))]     %A' 
Cprime = [0.007, -1]                  %VOC(SOC);
Eprime = [1, 0; 0, 1]              %E'
Fprime = [1, 0; 0, 1]              %F'

%
% Coefficients for probability
%
Rk = 1E-4;
RkP = 1E-4;
Qk1 = [2.5*10^(-7), 0; 0, 0];

%
% Initializing xhat and P (covariance matrix)
%
arrayOfxhats = zeros(2, totalTime);   
rawxhats = zeros(2, totalTime);
arrayOfxhats(1) = 1;                        %Assumes SOC starts at 1 and Vc starts at 0
LKarray = zeros(2, totalTime); 
diff = zeros(2, totalTime);
arrayOfPs = zeros(2, 2, totalTime);
arrayOfPs(1:2, 1:2, 1) = RkP * eye(2);       %P starts of as n 2x2 identity matricies
KGP = zeros(1, totalTime);

%--------------------------------------------------------------------------
% At time k
%--------------------------------------------------------------------------

%Variables Qualitatively                          
%     (Constant) VOC: Open circuit voltage
% (Time variant) Vc: Voltage across capacitor in battery circuit model
% (Time variant) V: Voltage we measure (output from battery circuit model)
%     (Constant) Ccap: Capacitor's capacitance
%     (Constant) Rc: Resistor in parallel with capacitor
%     (Constant) R0: Ohmic resistance
%     (Constant) Cbat: Capacity of battery in AmpHours?
% 

%
% Function we need to calculate
%
VOC = @(SOC, Voc0) 0.007*SOC + Voc0; %????               %Need equation for voc in terms of SOC 

%yk = @(V, Voc0) V - Voc0;                         %yk: Actual Voltage
yk = @(V, Vc) V + Vc;   

hk = @(SOC, I, Voc0) VOC(SOC*100, Voc0) - R0*I;           %hk: Cacluated Voltage?

fk = @(xhatk_1, I, dt, Cbat, Ccap, Rc) xhatk_1 + dt * [-I / Cbat; (I / Ccap)  - (xhatk_1(2, 1) / (Ccap * Rc))];  %xhat = previous xhat + change in xhat in time: dt     %XHAT IS A DERIVATIVEEEE

%
%Iterate for each sample
% t: Time
% totalTime
for t = 1:totalTime-1                                            %Not sure if this I is offset
    [arrayOfxhats(:, t+1), arrayOfPs(:, :, t+1), LKarray(:, t), diff(:, t), rawxhats(:, t), KGP(t)] = EKF(arrayOfxhats(:, t), arrayOfPs(:, :, t), I(t+1), I(t), V(t+1), Voc0, Rk, Aprime, Cprime, Eprime, Fprime, fk, dt, Cbat, Ccap, Rc, Qk1, yk, hk, actualSOC(t));

end

%Plotting
figure
%plot(timeSteps, actualSOC, 'DisplayName', "Actual SOC")
allsteps = (V + I*R0 + Vc - Voc0)/0.7;
allsteps(1) = 1;
plot(timeSteps, allsteps , 'DisplayName', "Actual SOC")
hold on
plot(timeSteps, arrayOfxhats(1, :), 'DisplayName', "Extended Kalman Filter")
%plot(timeSteps, (V + I*R0 + Vc - Voc0)/0.007, 'DisplayName', "Extended Kalman Filter")
%plot(t, SOCdr, 'DisplayName', "Open Loop")
ylim = ([72-20, 72+20]);
xlim([0, totalTime])
legend
title(["SOC Estimate Using Extended Kalman Filter, Open Loop Estimate and", "Actual SOC vs. Time"])

% Page 14
xlabel("time (s)")
ylabel("SOC")
saveas(gcf, "./Figures/2ekf.jpg")

%Vc should change with time??


%Calculate next xhat and P using the previous ones
function [xhatCorrected, PCorrected, Lknum, diff, hy, KGP] = EKF(xhatk_1, Pk_1, I, Ik_1 , V, Voc0, Rk, Aprime, Cprime, Eprime, Fprime, fk, dt, Cbat, Ccap, Rc, Qk1, yk, hk, actualSOC)

    %Testing
    %xhat gives us current SOC
    hy = xhatk_1; % - dt * [I / Cbat; 0];
    %End Testing


    xhat = fk(xhatk_1, I, dt, Cbat, Ccap, Rc);
    P = Aprime * Pk_1 * Aprime.' + Eprime * Qk1 * Eprime.';

    %Testing
    KGP = Cprime * P * Cprime.' + Rk;
    Lknum = P * Cprime.';
    %End Testing

    Lk = P * Cprime.' * (Cprime * P * Cprime.' + Rk)^-1;
    %                            Measured Voltage       What the voltage
                                                      % Should be at SOC
    
    %Testing                                               
    diff = yk(V, xhat(2, 1)) - hk(xhat(1, 1), I, Voc0);
    %xhatCorrected = xhat + Lk * (yk(V, Voc0) - hk(xhat, I, Voc0));
    xhatCorrected = xhat + Lk * (yk(V, xhat(2, 1)) - hk(xhat(1, 1), I, Voc0));
    %End Testing


    %xhatCorrected = xhat + Lk * (hk(actualSOC, I, Voc0) - hk(xhat(1, 1), I, Voc0));
    PCorrected = P - Lk * Cprime * P;
end


