% This file contains a brief description of all the files that are included
% in this package. A more detailed description is available in each
% individual file. All files are intended for use with Matlab.
%
% Data sets:
%    FurnaceHeating.mat   - Temperature measurements 5, 11 and 17 mm below
%                           the surface of a steel slab. The data are
%                           recorded at 10 Hz and collected in three vectors
%                           TC1, TC2, and TC3.
%   
%   ThermalConductivity.m - The thermal conductivity k(T) for the specific 
%                           steel composition in the experiment.
%   
%   Density_SpecificHeat.m - The product of the density and the specific
%                            heat capacity for the steel composition in the
%                            experiment. 
%
%
% Main functions:
%
%   SHEForwardSolv.m    - Finite difference solver for the heat equation with 
%                         temperature dependent material properties. This
%                         is used to create numerical test problems as well
%                         as computing thermal gradients from temperature 
%                         measurements. This code solves a well-posed
%                         problem.
%  
%  SHESolver.m - Solves the Cauchy problem for the heat equation. This is
%                an ill-posed problems. Regularization is done either by
%                the Fourier transform or by Cubic smoothing splines. 
%
%  PaperGraphs.m - This file runs all the experiments presented in our 
%                  paper "Estimation of the Surface Heat Flux using Cubic
%                  Smoothing Splines. The print commands are commented out
%                  so no postscript files are produced.
%
% Utility routines                 
%  
%  SHESystem.m - The inverse problem is solved by rewriting it as a system
%                of ODEs. This function contains the right-hand-side function 
%                for the ODE. It is called by SHESolver.m
%
%  SSPDeriv.m - Computation of the time derivative by using Cubic smoothing 
%               splines. 
% 
%  FFTDeriv.m - Computation of the time derivative by the Fourier
%               transform. This function is not used by included for
%               completeness.
%
%
% By: Fredrik Berntsson, Linköping university, Sweden. 
%     Email: fredrik.berntsson@liu.se
