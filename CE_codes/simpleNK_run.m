% simple examples
clear all
clc
close all

%
syms c pi cp pip ups eps upsp epsp
syms sig_1 alfa beta kapa rhoe rhou

f0     = - pi + beta*pip + kapa*c + ups;
f1(1,1)= - c + cp - sig_1*(alfa*pi - pip) + eps;
f1(2,1)= - upsp + rhou*ups;
ftheta = - epsp + rhoe*eps;

f = [f0; f1; ftheta];

x      = ups;
xp     = upsp;
theta  = eps;
thetap = epsp;
y      = c;
yp     = cp;
z      = pi;
zp     = pip;

control  = [y  z ];
controlp = [yp zp];
state    = [x  theta];
statep   = [xp thetap];

ntheta = length(theta);
nTHETA = sum(ntheta);
nz = length(z);
nZ = sum(nz);
nx = length(x);
ny = length(y);

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
GX
nRHO =  -inv(nfxp(end-nTHETA + 1:end,nx + 1:nx + nTHETA))*nfx(end-nTHETA + 1:end,nx + 1:nx + nTHETA);%
[ hX, gX ] = gx_hx_partinfo(nRHO, GX, HX, nfx, nfy, nfxp, nfyp,  nz, ntheta);
gX
% Doublecheck the formulas in the text
D = - kapa*sig_1*(alfa-rhoe)- (rhoe-1)*(beta*rhoe-1);
check_GX = [ (alfa-rhoe)*sig_1 beta*rhoe-1; rhoe-1 -kapa ]/D
check_gX = [ (alfa-rhoe)*sig_1  beta*rhoe-1-alfa*sig_1*kapa      alfa*sig_1*rhoe*kapa;...
            rhoe-1   0   -kapa*rhoe ]/D