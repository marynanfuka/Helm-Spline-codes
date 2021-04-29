% SHEForwardSolv: Solution to the direct heat equation problem, 
%   that is given boundary data T(x=L0,t)=f(t) we compute the
%   derivative Tx along the line x=L0 and x=L1. The equation that is 
%   solved is
%
%      (a(T)*Tx)x=b(T)*Tt,  L0<x<L1, 0<t<Tend.
%      T(L0,t)=F(t),     0<t<Tend.
%      T(L1,t)=G(t),     0<t<Tend.
%
% Usage:
%
% >>[H,W]=SHEForwardSolv(x,t,F,G,a,b); 
%
% Input parameters are
%  x,t - The grids in time and space direction respectively. Its assumed
%        the known surface is x(1)=L0.
%  F   - known temperature T(L0,t)=F.
%  G   - known temperature T(L1,t)=G.
%  a,b - Coefficients in the equation.
%  
% Output parameters:
%
%  H,W - The derivatives Tx(L0,t) and Tx(L1,t). 
%
function [H,W,Tsurf]=SHEForwardSolv(x,t,F,G,a_fun,b_fun)

%
% Compute size of time and space grids and also the step lengths.
%
 [n1,n2]=size(x);,if n1<n2,x=x';,end
 [n1,n2]=size(x);,if n1<n2,t=t';,end
 
 N=length(t);dt=t(2)-t(1);
 M=length(x);dx=x(2)-x(1);
 
%
% Check if a_fun and b_fun are in fact scalars. 
%
 if isnumeric(a_fun),a_fun=@(T)a_fun*ones(size(T));,end
 if isnumeric(b_fun),b_fun=@(T)b_fun*ones(size(T));,end

 
%
% Solve the stationary problem for t=0. 
%
 [T]=StationaryTemperature(x,a_fun,F(1),G(1));
 
%
% Create vectors to store the output vectors. Also compute
% the output at t=0.
%
 H=zeros(size(G));W=H;
 H(1)=(-3*T(1)+4*T(2)-T(3))/2/dx; 
 W(1)=(+3*T(M)-4*T(M-1)+T(M-2))/2/dx;
 
if nargout>2,Tsurf=zeros(N,M);Tsurf(1,:)=T';,end
%
% Perform timesteps
%
 for j=2:N
     
     [T]=OneTimeStep(x,dt,T,F(j),G(j),a_fun,b_fun  );
     
     % Comput output
     H(j)=(-3*T(1)+4*T(2)-T(3))/2/dx; 
     W(j)=(+3*T(M)-4*T(M-1)+T(M-2))/2/dx;
     if nargout>2,Tsurf(j,:)=T';,end
 end

 
end

%======================================================================%
%
% Solve the stationary problem (a(T)Tx)x=0 using fixed point iterations.
%
function [T]=StationaryTemperature(x,a_fun,F,G)

 % Initial temperature using linear interpolation
 p1=polyfit([x(1) x(end)],[F G],1);T=polyval(p1,x);
 m=length(x);
 
 % Fixed-point iteration: Evaluate a(T) -> Solve (a(x)Tx)x=0 -> repeat 
 for j=1:5
   a=a_fun( (T(2:end)+T(1:end-1))/2 );
   A=spdiags([a(1:m-2) -(a(1:m-2)+a(2:m-1)) a(2:m-1)],0:2,    m-2,m);
   T=full([F;  A(:,2:m-1)\(-F*A(:,1)-G*A(:,m)) ; G]);
 end
 
end


%======================================================================%
%
% Local subroutine to perform one time-step with stepsize dt.
%
function [T_next]=OneTimeStep(x,dt,T,F_next,G_next,a_fun,b_fun  )


 % Create the discretization of the space derivatives (aTx)x. 
 m=length(x);dx=x(2)-x(1);
 a=a_fun( (T(2:end)+T(1:end-1))/2 );b=b_fun(T(2:m-1));
 A=spdiags([a(1:m-2) -(a(1:m-2)+a(2:m-1)) a(2:m-1)],0:2,m-2,m)*(dt/2/dx^2);
 I=spdiags(b,0,m-2,m-2);
 T_next=(I-A(:,2:m-1))\(I*T(2:m-1)+A*T+F_next*A(:,1)+G_next*A(:,m));
 T_next=[F_next;T_next;G_next];
 
end
