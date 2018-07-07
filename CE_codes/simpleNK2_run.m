% simple examples
clear all
clc
close all

%
syms c pi cp pip ups eps upsp epsp
syms sig_1 alfa beta kapa rhoe rhou

f0(1,1)     = - pi + beta*pip + kapa*c + ups;
f1(1,1 )    = - c + cp - sig_1*(alfa*pi - pip) + eps;
ftheta(1,1) = - upsp + rhou*ups;
ftheta(2,1) = - epsp + rhoe*eps;

f = [f0; f1; ftheta];

theta1 = ups; theta1p = upsp;
theta2 = eps; theta2p = epsp;
theta  = [theta1 theta2];
thetap = [theta1p theta2p];

z0      = pi ; z0p     = pip;
z1      = c ;  z1p     = cp;

z        = [z1  z0];
zp       = [z1p z0p];
control  = [ z ];
controlp = [ zp];
state    = [ theta];
statep   = [ thetap];

ntheta = [length(theta1); length(theta2) ];
nTHETA = sum(ntheta);
nz0    = length(z0);
nz1    = length(z1);
nz     = [nz1; nz0];
nZ  = sum(nz);
nx  = 0;
ny  = 0;

approx = 2;
[ fx, fxp, fy, fyp, fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,state,control,statep,controlp,approx);

ups   = 0; eps = 0; pi = 0; c = 0;
upsp  = ups; epsp = eps; pip = pi; cp = c;
sig_1 = 1;
alfa  = 1.5;
beta  = 0.99;
kapa  = 0.1;
rhoe  = 0.9; 
rhou  = 0.9; 

num_eval()

[GX,HX,outcome] = gx_hx(nfy, nfx, nfyp, nfxp);
nRHO =  -inv(nfxp(end-nTHETA + 1:end,nx + 1:nx + nTHETA))*nfx(end-nTHETA + 1:end,nx + 1:nx + nTHETA);
[ hX, gX ] = gx_hx_partinfo(nRHO, GX, HX, nfx, nfy, nfxp, nfyp,  nz, ntheta);
gX
% Doublecheck formulas in the text
D = -kapa*sig_1*( alfa - rhoe)- (rhoe-1)*(beta*rhoe-1);%(alfa-rhoe)*sig_1+(rhoe-1)*alfa*sig_1
check_gX = [sig_1*(alfa-1)*rhoe 0  -alfa*sig_1*rhoe*(rhoe-1)  rhoe*(beta*rhoe-1) ; ...
    0 0 (rhoe-1)*rhoe -rhoe*kapa]/D