%SHESolver: Solve the sideways heat equation, 
%
% ( kappa Tx )x = Tt,                 L1 < x < L2,
% T(L2,t)=G(t), Tx(L2,t)=H(t),        0 < t < Tend,
% Tt(x,0)=0,                          L1 < x < L2.
%
% Where the initial condition means the material is initally at steady 
% state.
%
% Usage:
%  >>  [T,Tx,x]=SHESolver( x , t , G , H ,  kappa , regu_method , regu_param );
%
% where 
%  x - is a vector with the space coordinates where the solution is sought.
%      L1=x(1) and L2=x(end). 
%  t - is a vector containing the time grid such that t(1)=0 and
%      t(end)=Tend. 
%  G,H - and vectors containing the values g(t) and h(t) respectively.
%
%  regu_method - Method used to approximate the time derivative. Is either
%     'SSP' or 'FFT'.
%  regu_param - regularization parameter for the medtod used.
%
function [U]=HelmholtzCauchy(x , y , F1 , H1 , k2 ,regu_method,regu_param);

%
% Make sure G and H are column vectors
% 
 [n1,n2]=size(F1);if n1<n2,F1=F1';,end;
 [n1,n2]=size(H1);if n1<n2,H1=H1';,end;
 F1(1)=0;F1(end)=0;H1(1)=0;H1(end)=0;

% Initial filtering of F1,H1.
%
 if strcmp(regu_method,'SSP')
  F1=SSPDeriv(F1,regu_param,0); 
  H1=SSPDeriv(H1,regu_param,0);
 else
  error('HelmholtzCauchy: Unknown regularization method');
 end
 F1(1)=0;F1(end)=0;H1(1)=0;H1(end)=0; 
 ode_params=odeset('NormControl','on','RelTol',1e-4);
 [y,Z]=ode23('HelmholtzSystem', y , [F1;H1] ,ode_params,k2,regu_method,regu_param);
% 
% Collect output
%

  U=Z(:,1:length(F1));

