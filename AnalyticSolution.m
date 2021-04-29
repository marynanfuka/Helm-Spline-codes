function [F1,H1,F2,H2,U]=AnalyticSolution( x , y , L , k2 ); 

     
% Create the Dirichlet data u(x,L)=F2(x).

f_fun=@(x,y) sin(pi*x).*cosh(sqrt(pi^2-k2)*y)+...
             sin(2*pi*x).*sinh(sqrt(4*pi^2-k2)*y)+...
             sin(3*pi*x).*cosh(sqrt(9*pi^2-k2)*y); 
h_fun=@(x,y) sqrt(pi^2-k2)*sin(pi*x).*sinh(sqrt(pi^2-k2)*y)+... 
             sqrt(4*pi^2-k2)*sin(2*pi*x).*cosh(sqrt(4*pi^2-k2)*y)+...
             sqrt(9*pi^2-k2)*sin(3*pi*x).*sinh(sqrt(9*pi^2-k2)*y); 

F1=f_fun(x,0);F2=f_fun(x,L);H1=h_fun(x,0);H2=h_fun(x,L);


if nargout>4
  [X,Y]=meshgrid(x,y);   
  U =  sin(pi*X).*cosh(sqrt(pi^2-k2)*Y)+... 
       sin(2*pi*X).*sinh(sqrt(4*pi^2-k2)*Y)+...
       sin(3*pi*X).*cosh(sqrt(9*pi^2-k2)*Y); 
end