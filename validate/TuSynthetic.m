function [U, t, x] = TuSynthetic( varargin )
%function Tu2014
%TU2014 Synthetic dataset from the Tu et al (2014) Exp Fluids
%
% Tu, Jonathan H., Clarence W. Rowley, J. Nathan Kutz, and Jessica
% K. Shang. 2014. “Spectral Analysis of Fluid Flows Using Sub-Nyquist-Rate PIV
% Data.” Experiments in Fluids 55 (9): 1–13. doi:10.1007/s00348-014-1805-6.
%

% parameters
omega1 = 1.30;
omega2 = 8.48;
dt = 0.05;

% grid
x = linspace(-2,2,200);
y = linspace(-2,2,100);

% spatial mode fields
v1 = @(x,y) 2*exp( -(x-0.5).^2/(2*0.6^2) - (y-0.5).^2/(2*0.2^2) );
v2 = @(x,y) 2*exp( -(x+0.25).^2/(2*0.6^2) - (y-0.35).^2/(2*1.2^2) );

% evaluate spatial mode fields on a grid
V1 = bsxfun(v1, x, y');
V2 = bsxfun(v2, x, y');
[X,Y] = meshgrid(x,y);

subplot(1,2,1); pcolor(X,Y,V1); axis square; shading interp
subplot(1,2,2); pcolor(X,Y,V2);axis square; shading interp

V1 = V1(:);
V2 = V2(:);

t = 0:dt:400;

Osc1 = sin(omega1*t);
Osc2 = sin(omega2*t);

S = bsxfun( @times, V1, Osc1 ) + bsxfun( @times, V2, Osc2 );
S = bsxfun( @plus, S, rand(1,numel(t)) );
