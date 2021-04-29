
function [Kappa]=ThermalConductivity(T);

%  Table over values k(T);

  kvals= [65.3 62.8 60.3 57.8 55.7 53.2 51.1 48.6 46.5 43.5 41.0 39.4 37.7 36.0 33.9 31.8 30.1 27.6 27.2 27.2 27.6 28.1 28.5 29.3 29.7 30.9];
  tvals= [0 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1300];
 

% Compute interpolating spline 

  pp=csape(tvals,kvals);

% splinefunktionen är definierad för min(tvals)<T<max(tvals)

  tmp=kvals(1)*ones(size(T));
  ind=find( T > max(tvals) ); tmp(ind)=kvals(length(kvals));
  ind=find( min(tvals) <= T & T <= max(tvals) );
  tmp(ind)=ppval(pp,T(ind));
  Kappa=tmp;
 