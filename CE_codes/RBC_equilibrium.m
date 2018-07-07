function [y, z, x, theta, fx, fxp, fy, fyp, fypyp, fypy,  fypxp,  fypx,  fyyp,  fyy,  fyxp,  fyx,  fxpyp,  fxpy,  fxpxp,  fxpx,  fxyp,  fxy,  fxxp, fxx , f, nz, ntheta] = RBC_equilibrium;
% This program computes the equilibrium vector and its derivatives for the
% RBC model with labor and timing constraints.
% timing constraint: Technology shock realizes in the middle of the period.
% Capital and investment must be chosen before the shock realizes.

%(c) Anna Kormilitsina, June 2010
global approx
syms ALFA SIG PHI BETTA DELTA RHOa RHOksi
syms k_cu    a_cu   h_cu   c_cu   i_cu  ksi_cu
syms k_cup   a_cup  h_cup  c_cup  i_cup ksi_cup
syms k_ba1  
syms k_ba1p 
syms K       A      H      C      I KSI

ksi_cu = KSI;
ksi_cup = KSI;

prf       = exp(a_cu) * k_ba1^ALFA *h_cu^(1-ALFA);
prfp      = exp(a_cup)* k_cu^ALFA*h_cup^(1-ALFA);
dprf_dh   = diff(prf, h_cu);
dprfp_dkp = diff(prfp, k_cu);

util_cu   = ksi_cu *log(c_cu ) - h_cu^2/2;
util_cup  = ksi_cup*log(c_cup) - h_cup^2/2;

dup_dcp   = diff(util_cup, 'c_cup');
dudc      = diff(util_cu , 'c_cu' );
dudh      = diff(util_cu , 'h_cu' );

% Eq 1 - kp, partial info
f_0(1,1)  = - dudc + BETTA * dup_dcp * (dprfp_dkp + 1 - DELTA);
% Eq 2 - c, full info, because prf is full info since it depends on h
f_1(1,1)  = - c_cu + prf  - k_cu  + ( 1 - DELTA) * k_ba1;
% Eq 3 - h, full info
f_1(2,1)  = dudh  + dudc*dprf_dh;
f_1(3,1) = + k_cu - k_ba1p;
%==========================================================================
% second subperiod shock
f_theta1(1,1) = - (a_cup) + RHOa * (a_cu);
% Partial-info controls
z0      = [  k_cu  ];
z0p     = [  k_cup ];
Z0      = [  K     ];
nz0     = length(z0);

z   = [ z0  ];
zp  = [ z0p ];
Z   = [ Z0  ];
nz  = [ nz0 ];

% Full-info controls
y      = [c_cu   h_cu  ];
yp     = [c_cup  h_cup ];
Y      = [C      H     ];

ny     = length(y);

% Endogenous states + beginning period shocks
x      = [ k_ba1  ];
xp     = [ k_ba1p ];
nx     = length(x);
X      = [ K      ];
% Middle-period shocks

theta0      = [];
theta0p     = [];
THETA0      = [];
x           = [ x  ];
xp          = [ xp ];
X           = [ X  ];
nx          = length(x);

theta1      = [ a_cu  ];
theta1p     = [ a_cup ];
THETA1      = [ A     ];

theta       = [ theta1  ];
thetap      = [ theta1p ];
THETA       = [ THETA1  ];
ntheta      = [length(theta1) ]; % vector showing number of shocks with realization in a subperiod

nTHETA     = sum(ntheta);
nM         = length(ntheta);

%Create function f
f = [ f_0; f_1; f_theta1 ];
vs = [ x, theta, y, z, xp, thetap, yp, zp];

%f = subs(f, vs, exp(vs));
[ fx, fxp, fy, fyp, fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,[x , theta],[y, z],[xp thetap],[yp zp],approx); 
