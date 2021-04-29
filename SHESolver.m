%SHESolver: Solve the sideways heat equation, 
%
% ( a(T) Tx )x = b(T)*Tt,             L1 < x < L2,
% T(L2,t)=G(t), Tx(L2,t)=H(t),        0 < t < Tend,
% Tt(x,0)=0,                          L1 < x < L2.
%
% Where the initial condition means the material is initally at steady 
% state.
%
% Usage:
%  >>  [T,Tx,x]=SHESolver( x , t , G , H ,  a_fun , b_fun , regu_method , regu_param );
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
function [T,Tx,x]=SHESolver( x , t , G , H , a_fun , b_fun , regu_method , regu_param );


%
% Make sure G and H are column vectors
% 
 [n1,n2]=size(G);if n1<n2,G=G';,end;
 [n1,n2]=size(H);if n1<n2,H=H';,end;

%
% Check if a_fun and b_fun are in fact scalars. 
%
 if isnumeric(a_fun),a_fun=@(T)a_fun*ones(size(T));,end
 if isnumeric(b_fun),b_fun=@(T)b_fun*ones(size(T));,end

%
% Redefine a_fun to include Tend from the change of variables to 0<t<1.
%
 Tend=t(end)-t(1);
 a_fun2=@(T)Tend*a_fun(T);

%
% Initial smoothing
%

%
% Wrong: Have to adjust regularization paramater when [0,Tend]->[0,1]
%
 if strcmp(regu_method,'FFT')
   regu_param=regu_param*Tend;
   [~,G]=FFTDeriv(G,regu_param);
   [~,H]=FFTDeriv(H,regu_param);
 elseif strcmp(regu_method,'SSP')
   [~,G]=SSPDeriv(G,regu_param);
   [~,H]=SSPDeriv(H,regu_param);
 else
    error('SHESolver: Unknown regularization method');
 end
 

%
% Solve problem using ode45. Parameters kappa and xi_c are passed to the
% ODE system function.
%
 [x,Z]=ode45('SHESystem', x(end:-1:1) , [G;a_fun2(G).*H] , [] , a_fun2 , b_fun , regu_method , regu_param );
 
% 
% Collect output
%
 T=Z(:,1:length(G));Tx=Z(:,length(G)+1:end)./a_fun2(T);
 
