function sol = ComplexNLSQ( f, guess )
%COMPLEXOPTIMIZE
% Minimize a (vecor) function f of complex variables using Matlab's
% standard (real) optimizers.

  CtoR = @(z)[real(z); imag(z)];
  RtoC = @(x)complex(x(1:(numel(x)/2)),x((numel(x)/2)+1:end));

  tt = complex( rand(10,1), rand(10,1) );

  assert( all( RtoC(CtoR(tt)) == tt ) )

  fR = @(x) f( RtoC(x) );

  opts = optimset(@lsqnonlin);
  opts.Display = 'iter';

  [x,resnorm,residual,exitflag,output] = lsqnonlin(fR,CtoR(guess),[],[],opts);

  sol = RtoC(x);

end
