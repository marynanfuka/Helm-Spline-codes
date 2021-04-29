function [RhoCp]=Density_SpecificHeat(T);

%  Tabell med värden rho(T)*Cp(T)
   
 A=10^5;
 rho_cp= A.*[35.7 37.8 39.3 40.6 41.8 43.1 44.6 46.1 48.3 51.1 54.5 57.9 61.2 66.2 84.2 66.5 60.4 64.3 50.3 50.7 50.5 50.3 50.1 49.9 49.8 49.6 49.5];
 tvals= [0 75 125 175 225 275 325 375 425 475 525 575 625 675 725 775 825 875 925 975 1025 1075 1125 1175 1225 1275 1300];


% interpolera med en naturlig spline  

 pp=csape(tvals,rho_cp,'variational');

% splinefunktionen är definierad för min(tvals)<T<max(tvals)

  tmp=rho_cp(1)*ones(size(T));
  ind=find( T > max(tvals) ); tmp(ind)=rho_cp(length(rho_cp));
  ind=find( min(tvals) <= T & T <= max(tvals) );
  tmp(ind)=ppval(pp,T(ind));
  RhoCp=tmp;

