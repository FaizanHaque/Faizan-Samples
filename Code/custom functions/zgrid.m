function [zt,zw,dzt,dzw] = zgrid(nz);

%function [zt,zw,dzt,dzw] = zgrid(nz,H);
% [zt,zw,dzt,dzw] = czgrid(km);
% MATLAB script to create a smooth vertical grid for the OGCM
% nz wet points with 
%    one extra dry air point at the surface and 
%    one extra dry solid earth point at the bottom
% Author: Francois Primeau 13/10/2000
% 
%
% imput parameter------------------------------------------------------
%

  s = 0.3609;     %stretching factor
  H = 5750;  %maximum depth 
  formula = 'tanh((1-eta)/s)/tanh(1/s)-1'; %smooth transformation
  result = '2*(2*cosh(((eta-1)^2)/s)-3)/(s^2*cosh((eta-1)/s)^2)';
  % ---------------------------------------------------------------------
  
  %
  % compute the zt,zw,dzt,dzw and metric factor in truncation error------
  %
  km = nz;
  eta = -1/km:1/km:1+1/km;
  f = inline(vectorize(formula),'eta','s');
  ef = inline(vectorize(result),'eta','s');
  zw = H*(f(eta,s));
  
  dzt = zw(2:end) - zw(1:end-1);
  dzw = 0.5*(dzt(2:end)+dzt(1:end-1));
  zt = (zw(1)+dzt(1)/2)+cumsum([0,dzw(1:length(dzw))]);

  zt = zt(1:end);
  zw = zw(1:end-1);
  dzt = abs(dzt(1:end));
  dzw = abs([dzw,dzw(end)]); % the bottom grid cell is always land so its
                        % thikness is irrelevant



