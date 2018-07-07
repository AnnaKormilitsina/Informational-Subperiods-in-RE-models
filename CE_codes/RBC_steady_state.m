function [DELTA,ALFA,BETTA,RHOa, RHOksi, A, KSI, H,K,C,I] = RBC_steady_state
% This program computes the steady state for the
% RBC model with labor.

%(c) Anna Kormilitsina, June 2010

BETTA  = 1.03^(-0.25); %discount rate
DELTA  = 0.025; %depreciation rate
ALFA   = 0.3; %capital share
%PHI    = 0.7; % utility parameter
RHOa   = 0.9; %persistence of technology shock
RHOksi = 0.9;
%SIG    = 0; %intertemporal elasticity of substitution

KSI    = 1;
A      = 0; %steady-state value of technology shock 
K_H    = ((1/BETTA+DELTA-1)/exp(A)/ALFA)^(1/(ALFA-1)); %steady-state value of capital
I_H    = DELTA *K_H;
C_H    = exp(A)*K_H^ALFA-DELTA*K_H; 
%H     = (1 + C_H/(1/PHI-1)/(1-ALFA)/A/K_H^ALFA )^(-1) ;
H      = (exp(A)*K_H^ALFA*KSI*(1-ALFA)/C_H)^(1/2);
K      = K_H*H;
I      = I_H*H;
C      = C_H*H;