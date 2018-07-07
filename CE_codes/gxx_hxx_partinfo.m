function [gXX,hXX]  = gxx_hxx_partinfo(nfX,nfXP,nfY,nfYP,nfYPYP,nfYPY,nfYPXP,nfYPX,nfYYP,nfYY,nfYXP,nfYX,nfXPYP,nfXPY,nfXPXP,nfXPX,nfXYP,nfXY,nfXXP,nfXX,hX,gX,GXX,HXX, nz, ntheta, P);
% This code can be used to obtain the second-order dynamics of the 
% RE model with timing constraints
% from the dynamics of the full information version of the model.
% The algorithm is described in "Solving Rational Expectations Models 
% with Partial Information Structure: A Perturbations Approach" by
% Anna Kormilitsina

% Outputs:
% hXX - 3-dimensional array of size nx by (nx + 2*sum{ntheta)) by (nx + 2*sum{ntheta))
% gXX - 3-dimensional array of size (nY+nZ) by (nx + 2*(ntheta) ) by (nx + 2*(ntheta) )
% hXX and gXX provide the second-order approximations to the dynamics of the model 
% state and control vectors in the following form:
% [xp; thetap ] = hX * [ x; theta; theta(-1)] +  0.5*hXX *_2[ x; theta; theta(-1) ]^T *_3 [ x; theta; theta(-1) ]^T +0.5*hSS*sigma^2 + sigma epsilon'
% and
% [y; z] = gX * [ x; theta; theta(-1) ] + 0.5*gXX *_2[ x; theta; theta(-1) ]^T *_3 [ x; theta; theta(-1) ]^T  + 0.5 gSS*sigma^2,
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
% HX,gX,GXX,HXX, nz, ntheta, eqs, P
% hX and gX are the first-order approximations to the dynamics of the model 
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
% P -  is the variance-covariance matrix of shocks with realization in the
% middle of the period (must be block-diagonal, with blocks corresponding 
% to information sets),
% nfX,nfXP,nfY,nfYP,nfYPYP,nfYPY,nfYPXP,nfYPX,nfYYP,nfYY,nfYXP,nfYX,
%   nfXPYP,nfXPY,nfXPXP,nfXPX,nfXYP,nfXY,nfXXP,nfXX, - first and second 
%     order derivatives of a vector of equilibrium conditions f,
%       evaluated at the steady state
% GXX,HXX - solution of the full-information version of the model
% nz - vector, each element i of which contains the number of control 
%      variables with decisions in information set i, 
% ntheta - vector, each element i of which contains  the number of shocks
%          with realizations in the beginning of information set i.
% eqs - indicator matrix; the number of rows equals the number of
%       subperiods-1, and the number of columns is equal to the number of
%       equations in the vector of equilibrium conditions f;
%       eqs(m,i) = 1 if equation f(i) is in the information set m, which
%       means this equation pins down some control variable in subperiod m.

% (c) Anna Kormilitsina, January 2011

nTHETA = sum(ntheta);
nZ     = sum(nz);               % number of all partial controls
nx     = size(nfX,2) - nTHETA;  % number of endogenous and exogenous states with realizations in the beginning of a period
ny     = size(nfY,2) - nZ;      % number of full info controls only
n      = nx + ny + nZ + nTHETA; % number of all variables and number of all equations
nM     = length(nz);            % number of subperiods

hx     = hX(1:nx, 1 : nx); 
Gx     = gX( :  , 1 : nx);
var    = [ Gx; eye(nx) ];
varx   = cat(1,   GXX(:,1:nx,1:nx),     zeros(nx,nx,nx));

hXX    = zeros(nx + 2*nTHETA,nx + 2*nTHETA,nx + 2*nTHETA);
gXX    = zeros(ny + sum(nz), nx + 2*nTHETA,nx + 2*nTHETA);    

htheta_1 = hX(1:nx,nx+nTHETA+1:nx+2*nTHETA);
varth_1  = [htheta_1; gX(:,nx+nTHETA+1:nx+2*nTHETA)];

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

   nfzm    = nfY(nzm+1:end-nTHETA,nym+1:ny+nZ);
   nfyzxp  = [nfYP(nzm+1:end-nTHETA,1:ny+nZ) nfXP(nzm+1:end-nTHETA,1:nx)];

   A       = [ nfyzxp*var, nfY(nzm+1:end-nTHETA,1:nym) ];
   invA    = inv(A);

   barjxthetam      = GXX(nym + 1:ny+nZ,1:nx, nx + nthetam_1 +1:nx + nthetam);
   jxthetam_1       = nmodeproduct3( barjxthetam , Pm',3 );
   jthetam_1x       = permute(jxthetam_1,[1 3 2]);
   barjthetamthetam = GXX(nym + 1:ny + nZ, nx+nthetam_1 + 1:nx + nthetam, nx + nthetam_1 + 1:nx + nthetam);
   jthetam_1thetam_1  = nmodeproduct3( nmodeproduct3( barjthetamthetam, Pm', 2), Pm' ,3);

   % step 2
   s2 = [ [ gX(:,1:nx+nTHETA+ntheta(m)); eye(nx+nTHETA+ntheta(m)) ]*hX(1:nx+nTHETA+ntheta(m),nx+1:nx+ntheta(m)); gX(1:nym,nx+1:nx+ntheta(m))]';
   nfYPy        = nfYPY(nzm+1:end-nTHETA,:,1:nym);
   nfYPtheta    = nfYPX(nzm+1:end-nTHETA,:,nx + nthetam_1 + 1: nx + nthetam); 
   dnfYPtheta   = nmodeproduct3( cat(3,nfYPYP(nzm+1:end-nTHETA,:,:), nfYPXP(nzm+1:end-nTHETA,:,:), nfYPtheta , nfYPy), s2, 3);% OK
   
   nfXPy        = nfXPY(nzm+1:end-nTHETA,:,1:nym);
   nfXPtheta    = nfXPX(nzm+1:end-nTHETA,:,nx + nthetam_1 + 1:nx + nthetam );
   dnfXPtheta   = nmodeproduct3( cat(3,nfXPYP(nzm+1:end-nTHETA,:,:), nfXPXP(nzm+1:end-nTHETA,:,:), nfXPtheta, nfXPy), s2 ,3); % OK
   dnfxptheta   = dnfXPtheta(:,1:nx,:);
   
   dnfyzxpthetam = cat( 2, dnfYPtheta, dnfxptheta  ); % OK
   prod         = nmodeproduct3(dnfyzxpthetam,  var(:,1:nx)', 2) + ...
        nmodeproduct3(nmodeproduct3( GXX(:,1:nx, 1: nx+nTHETA),hX(1:nx+nTHETA,nx + 1:nx+ntheta(m))',3), nfYP(nzm+1:end-nTHETA,1:ny+nZ) , 1 );
   
   nfYy         = nfYY(nzm+1:end-nTHETA,:,1:nym);
   nfYtheta     = nfYX(nzm+1:end-nTHETA,:,nx + nthetam_1 + 1: nx + nthetam);
   dnfYtheta    = nmodeproduct3( cat(3,nfYYP(nzm+1:end-nTHETA,:,:), nfYXP(nzm+1:end-nTHETA,:,:), nfYtheta, nfYy ), s2, 3); 
        
   B            = cat( 2, prod  , dnfYtheta);
   C            = - nmodeproduct3(  nmodeproduct3( B, varth_1(:,1:ntheta(m))', 2), invA ,1);

   htheta_1theta = C(1:nx,:,:);
   gtheta_1theta = C(nx+1:nx + nym,:,:);
   hthetatheta_1 = permute(htheta_1theta,[1 3 2]);
   gthetatheta_1 = permute(gtheta_1theta,[1 3 2]);

   %%% step 3
   s3 = [ Gx*htheta_1(:,1:ntheta(m)); varth_1(:,1:ntheta(m)) ]';
   dnfYPtheta_1   = nmodeproduct3( cat(3,nfYPYP, nfYPXP(:,:,1:nx), nfYPY), s3 , 3);
   dnfXPtheta_1   = nmodeproduct3( cat(3,nfXPYP, nfXPXP(:,:,1:nx), nfXPY), s3 , 3);
   dnfxptheta_1   = dnfXPtheta_1(:,  1:nx,  : );
   dnfyzxptheta_1 = cat( 2, dnfYPtheta_1, dnfxptheta_1);
   dnfyzxptheta_1 = dnfyzxptheta_1(nzm+1:end-nTHETA,:,:);
   
   nfYxp         = nfYXP(:,:,1:nx);
   dnfYtheta_1   = nmodeproduct3( cat(3,nfYYP, nfYxp, nfYY), s3 , 3); 
   dnfYtheta_1   = dnfYtheta_1(nzm+1:end-nTHETA,:,:);

   nfxYP = nfXYP(nzm+1:end-nTHETA, 1:nx,  :    );
   nfxY  = nfXY( nzm+1:end-nTHETA, 1:nx,  :    );
   nfxxp = nfXXP(nzm+1:end-nTHETA, 1:nx, 1:nx  );
   dnfxtheta_1   = nmodeproduct3( cat(3,nfxYP, nfxxp,nfxY), s3 , 3);
   dnfyzxtheta_1 = cat( 2, dnfYtheta_1, dnfxtheta_1);

   varx1nfyzxp = nmodeproduct3( varx , nfyzxp,1);
   
   %%%% Step 3 $h_{\theta_1,\theta_1}$ and $g_{\theta_1,\theta_1}$
   prod  = nmodeproduct3( dnfyzxptheta_1, var', 2) + ...
           nmodeproduct3( nmodeproduct3(GXX(:,1:nx,1:nx),  nfYP(nzm+1:end-nTHETA,:) ,1), htheta_1(:,1:ntheta(m))', 3 );
   prod1 = cat( 2, prod, dnfYtheta_1(:,:,1:ntheta(m)));
   B     = nmodeproduct3( prod1, varth_1(:,1:ntheta(m))' , 2); 
   C     = - nmodeproduct3( B , invA , 1 ) - nmodeproduct3( jthetam_1thetam_1, nfzm ,1);
    
   htheta_1theta_1 =  C(1:nx,:,:);
   gtheta_1theta_1 =  C(nx+1:nx + nym,:,:);

   % Step 4: hthetatheta and gthetatheta
   invPt = pinv(Pm'); 
   Hthetatheta = HXX(1:nx ,nx + nthetam_1 + 1:nx + nthetam, nx + nthetam_1 + 1 : nx + nthetam);
   Gthetatheta = GXX(1:nym,nx + nthetam_1 + 1:nx + nthetam, nx + nthetam_1 + 1 : nx + nthetam);
   hthetatheta = Hthetatheta - nmodeproduct3( hthetatheta_1  ,invPt,3) -...
       nmodeproduct3( permute(hthetatheta_1,[1 3 2]) , invPt,2) - nmodeproduct3(  nmodeproduct3( htheta_1theta_1, invPt,2)  , invPt,3);
   gthetatheta = Gthetatheta - nmodeproduct3( gthetatheta_1  ,invPt,3) - ...
       nmodeproduct3( permute(gthetatheta_1,[1 3 2]) , invPt,2) - nmodeproduct3(  nmodeproduct3( gtheta_1theta_1, invPt,2)  , invPt,3);

   % Step 5: gxtheta_1 and hxtheta_1
   prod  = nmodeproduct3( dnfyzxptheta_1(:,:,1:ntheta(m)), var', 2) + nmodeproduct3( varx1nfyzxp , htheta_1(:,1:ntheta(m))', 3 );
   prod1 = cat( 2, prod, dnfyzxtheta_1);
   B     = nmodeproduct3( prod1, [hx;var]', 2);
   C     = - nmodeproduct3( B , invA , 1 ) - nmodeproduct3( jxthetam_1(:,:,1:ntheta(m)), nfzm ,1);

   hxtheta_1 = C(1:nx,:,:);
   gxtheta_1 = C(nx+1:nx+nym,:,:);
   htheta_1x = permute(hxtheta_1,[1 3 2]);
   gtheta_1x = permute(gxtheta_1,[1 3 2]);
 
   % Step 6: gxtheta and hxtheta
   gxtheta = GXX(1:nym, 1:nx, nx+nthetam_1 + 1 : nx + nthetam) - nmodeproduct3(gxtheta_1, invPt, 3);
   hxtheta = HXX(1:nx , 1:nx, nx+nthetam_1 + 1 : nx + nthetam) - nmodeproduct3(hxtheta_1, invPt, 3);
   hthetax = permute(hxtheta,[1 3 2]);
   gthetax = permute(gxtheta,[1 3 2]);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   hXX(1:nx, 1:nx                             , nx + nthetam_1 + 1        : nx+nthetam                     ) = hxtheta;
   hXX(1:nx, 1:nx                             , nx + nthetam_1 + nTHETA+1 : nx+nTHETA+nthetam              ) = hxtheta_1;
   hXX(1:nx, nx + nthetam_1 + 1 : nx+nthetam   , nx + nthetam_1 + nTHETA+1 : nx+nTHETA+nthetam             ) = hthetatheta_1 ;
   hXX(1:nx, nx + nthetam_1 + nTHETA+1:nx+nTHETA+nthetam , nx +  nthetam_1 + nTHETA + 1 : nx+nTHETA+nthetam) = htheta_1theta_1 ;
   hXX(1:nx, nx + nthetam_1 + 1 : nx+nthetam             , 1            : nx                               ) = hthetax;
   hXX(1:nx, nx + nthetam_1 + 1 : nx+nthetam             , nx +  nthetam_1 + 1        : nx+nthetam         ) = hthetatheta ;
   hXX(1:nx, nx + nthetam_1 + nTHETA+1:nx+nTHETA+nthetam , 1            : nx                               ) = htheta_1x ;
   hXX(1:nx, nx + nthetam_1 + nTHETA+1:nx+nTHETA+nthetam , nx +  nthetam_1 + 1        : nx+nthetam         ) = htheta_1theta ;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   gXX(1:nym, 1:nx                                     , nx + nthetam_1 + 1        : nx+nthetam       )  = gxtheta;
   gXX(1:nym, 1:nx                                     , nx + nthetam_1 + nTHETA+1 : nx+nTHETA+nthetam)  = gxtheta_1;
   gXX(1:nym, nx+ nthetam_1+1:nx+nthetam               , 1                         : nx               )  = gthetax;
   gXX(1:nym, nx+ nthetam_1+1:nx+nthetam               , nx + nthetam_1 + 1        : nx+nthetam       )  = gthetatheta;
   gXX(1:nym, nx+ nthetam_1+1:nx+nthetam               , nx + nthetam_1 + nTHETA+1 : nx+nTHETA+nthetam)  = gthetatheta_1;
   gXX(1:nym, nx+ nthetam_1+nTHETA+1:nx+nTHETA+nthetam , 1                         : nx               )  = gtheta_1x;
   gXX(1:nym, nx+ nthetam_1+nTHETA+1:nx+nTHETA+nthetam , nx+ nthetam_1+1           : nx+nthetam       )  = gtheta_1theta ;
   gXX(1:nym, nx+ nthetam_1+nTHETA+1:nx+nTHETA+nthetam , nx+ nthetam_1+nTHETA+1    : nx+nTHETA+nthetam)  = gtheta_1theta_1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   gXX(nym+1:end, 1                     :nx               , nx+nthetam_1+nTHETA+1   : nx+nTHETA+nthetam) = jxthetam_1;
   gXX(nym+1:end, nx+ nthetam_1+nTHETA+1:nx+nTHETA+nthetam, 1                       : nx               ) = jthetam_1x;
   gXX(nym+1:end, nx+ nthetam_1+nTHETA+1:nx+nTHETA+nthetam, nx+ nthetam_1+nTHETA + 1: nx+nTHETA+nthetam) = jthetam_1thetam_1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   nym  = nym + nz(m);
   nzm  = nzm - nz(m);
end
gXX(1:ny+nZ          , 1:nx        , 1 : nx        ) = GXX(1:ny+nZ          , 1: nx, 1:nx);
hXX(1:nx             , 1:nx        , 1 : nx        ) = HXX(1:nx             , 1: nx, 1:nx);
hXX(nx+1:nx + nTHETA , 1:nx+nTHETA , 1 : nx+nTHETA ) = HXX(nx+1:nx + nTHETA ,  :   ,  :  );