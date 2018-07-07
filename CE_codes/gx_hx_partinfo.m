function [ hX, gX ] = gx_hx_partinfo(P, GX, HX, nfX,nfY, nfXP, nfYP, nz, ntheta )
% This code obtains the dynamics of the RE model with timing constraints
% from the dynamics of the full information version of the model.
% The algorithm is described in "Solving Rational Expectations Models 
% with Partial Information Structure: A Perturbations Approach" by
% Anna Kormilitsina

% Outputs:
% hX - square matrix of size nx + 2*sum{ntheta)
% gX - matrix of size ny+sum(nz) by nx + 2*(ntheta) 
% hX and gX are the first-order approximations to the dynamics of model 
% state and control vectors in the following form:
% [xp; thetap ] = hX * [ x; theta; theta(-1)] + epsilon'
% and
% [y; z] = gX * [ x; theta; theta(-1) ],
% where y, z, x, and theta are all in log-deviations from the steady state,
% and
%             x and xp - vectors of current and future state variables and 
%                      shocks with realization in the beginning of a period 
%                       (shocks with realization in the middle of the 
%                       period are not included in x)
%             y        - vector of full information control variables
%             z        - vector of partial information control variables
%             theta, thetap, and theta(-1) - vectors of current,
%                future, and previous period shocks with realization in the
%               middle of a period.
%              epsilonp - noise (not calculated)
% Inputs:
% P is the variance-covariance matrix of shocks with realization in the
% middle of the period (must be block-diagonal, with blocks corresponding 
% to information sets),
% GX and HX - matrices of linearized solution to the full-information
%             model
%  nfX,nfY, nfXP, nfYP,- first derivatives of a vector of 
%             equilibrium conditions f with respect to 
%             X, Y, XP, YP respectively, evaluated at the 
%             steady state
% nz - vector, each element i of which contains the number of control 
%      variables with decisions in information set i, 
% ntheta - vector, each element i of which contains  the number of shocks
%          with realizations in information set i.
% eqs - indicator matrix; the number of rows equals the number of
%       subperiods-1, and the number of columns is equal to the number of
%       equations in the vector of equilibrium conditions f;
%       eqs(m,i) = 1 if equation f(i) is in the information set m, which
%       means this equation pins down some control variable in subperiod m.

% (c) Anna Kormilitsina, January 2011

global nTHETA

nTHETA = sum(ntheta);
nZ     = sum(nz);
nx     = size(nfX,2)-nTHETA;
ny     = size(nfY,2)-nZ;
n      = nx + ny + nZ + nTHETA;
nM     = length(nz);
jx     = GX(ny+1:end,1:nx);
gx     = GX(1:ny, 1: nx);
hx     = HX(1:nx,1:nx);

nfy  = nfY(:,1:ny) ;
nfz  = nfY(:,ny+1:ny+nZ);
nfxp = nfXP(:,1:nx);

hX     = zeros(nx + 2*nTHETA,nx + 2*nTHETA);
gX     = zeros(ny + nZ, nx + 2*nTHETA);

nym = ny;
nzm = nZ;
for m = nM:-1:1
    nthetam = sum( ntheta(1:m));
    if m ~=1
        nthetam_1 = sum(ntheta(1:m-1));
    elseif m == 1
        nthetam_1 = 0;
    end
    Pm = P( nthetam_1 + 1 : nthetam, nthetam_1 + 1 : nthetam );
    jtheta_ba1 = GX( nym + 1:end, nx + nthetam_1 + 1 : nx + nthetam)*Pm; % jtheta_ba1^m
    gX( nym+1 : end, nx + nTHETA + nthetam_1 + 1 : nx + nTHETA + nthetam ) = jtheta_ba1;
    
    a = [ nfYP*GX(:,1:nx) + nfxp, nfY(:,1:nym) ];
    b = nfY(:,nym+1:end)*jtheta_ba1;
  
    a   = a(nzm+1:end-nTHETA, :);
    b   = b(nzm+1:end-nTHETA, :);
    hgtheta_ba1 = - inv(a)*b;
    htheta_1 = hgtheta_ba1(1:nx,:);
    gtheta_1 = hgtheta_ba1(nx+1:nx+nym,:);
    hX(1:nx,  nx + nTHETA + nthetam_1 + 1 : nx + nTHETA + nthetam) = htheta_1;
    gX(1:nym, nx + nTHETA + nthetam_1 + 1 : nx + nTHETA + nthetam) = gtheta_1;
    
    Gtheta   = GX(1:nym,nx + nthetam_1 + 1 : nx + nthetam);
    Htheta   = HX(1:nx ,nx + nthetam_1 + 1 : nx + nthetam);
    
    gX( 1:nym, nx + nthetam_1 + 1 : nx + nthetam ) = Gtheta - gtheta_1*pinv(Pm); % gtheta^m
    hX( 1:nx , nx + 1 : nx + nthetam)              = Htheta - htheta_1*pinv(Pm);% htheta^m
    % update nfy and nfz and the size of y vector
    nym                  = nym + nz(m);
    nzm                  = nzm - nz(m);
end
gX(1:ny,1:nx)     = gx;
gX(ny+1:end,1:nx) = jx;
hX(1:nx,1:nx)     = hx;

hX(nx+1       : nx+nTHETA  ,nx+1:nx+nTHETA) =  P;
hX(nx+nTHETA+1: nx+2*nTHETA,nx+1:nx+nTHETA) =  eye(nTHETA);