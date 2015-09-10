function sol = ComplexOptimize( f, guess )
%COMPLEXOPTIMIZE
% Minimize a (scalar) function f of complex variables using Matlab's
% standard (real) optimizers.

  CtoR = @(z)[real(z); imag(z)];
  RtoC = @(x)complex(x(1:(numel(x)/2)),x((numel(x)/2)+1:end));

  tt = complex( rand(10,1), rand(10,1) );

  assert( all( RtoC(CtoR(tt)) == tt ) )

  fR = @(x) f( RtoC(x) );

  opts = optimset(@fminsearch);
  opts.Display = 'iter';

  [x,fval,exitflag,output] = fminsearch(fR,CtoR(guess),opts);

  sol = RtoC(x);

end
