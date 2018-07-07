clear all
clc
close all

% This program computes impulse responses for an RBC model with labor and
% timing constraints using the algorithm described in "Solving Rational Expectations Models with Partial Information Structure: A Perturbations Approach"
% timing constraint: Technology shock realizes in the middle of the period.
% Capital and investment must be chosen before the shock realizes.

%(c) Anna Kormilitsina, January 2011

global approx
approx = 2;
T      = 5;

%Steady State and Parameter Values
[DELTA,ALFA,BETTA,RHOa, RHOksi, A, KSI, H,K,C,I] = RBC_steady_state;
%%%% solve the full-information version of the problem first
[y, z, x, theta, fx, fxp, fy, fyp, fypyp, fypy,  fypxp,  fypx,  fyyp,  fyy,  fyxp,  fyx,  fxpyp,  fxpy,  fxpxp,  fxpx,  fxyp,  fxy,  fxxp, fxx , f, nz, ntheta] = RBC_equilibrium;
var_cu      = [  y z x theta];
nh          = find( var_cu == 'h_cu' );
nc          = find( var_cu == 'c_cu' );
nk          = find( var_cu == 'k_ba1');
na          = find( var_cu == 'a_cu' );
%Evaluate derivatives of f under assumed calibration
h_cu   = (H); a_cu  = (A); k_cu = (K); c_cu = (C); i_cu = (I); ksi_cu = (KSI);
h_cup  = (H); a_cup = (A); k_cup = (K); c_cup = (C); i_cup = (I);ksi_cup = (KSI);
k_ba1p = (K); k_ba1 = (K); 

shock_size = 0.1;

num_eval()
nf

nx     = length(x);
ny     = length(y);
nz     = length(z);
ntheta = length(theta);
nTHETA = sum(ntheta);
nZ = sum(nz);

[GX,HX,outcome] = gx_hx(nfy, nfx, nfyp, nfxp);
[GXX,HXX] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,HX,GX);
eta   = zeros(nx+nTHETA,nTHETA);
state = [x theta];
eta(find(state == 'a_cu'),   find(theta == 'a_cu')) = shock_size ;
[Gss,Hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,HX,GX,GXX,eta);

if outcome == 0
    warning('problem here')
elseif outcome == 1
    %Now solve model with timing restrictions    
    nRHO =  -inv(nfxp(end-nTHETA + 1:end,nx + 1:nx + nTHETA))*nfx(end-nTHETA + 1:end,nx + 1:nx + nTHETA);
    [ hX, gX ] = gx_hx_partinfo(nRHO, GX, HX, nfx, nfy, nfxp, nfyp,  nz, ntheta);
    [gXX,hXX]  = gxx_hxx_partinfo(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hX,gX,GXX,HXX, nz, ntheta, nRHO);
end
GXX
gXX
HXX
hXX