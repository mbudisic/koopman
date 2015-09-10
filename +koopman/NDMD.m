function [lambda, M] = NDMD( snapshots, t, Nmd )
%NDMD An implementation of a Nonuniformly sampled DMD algorithm by Gueniat
%et al 2015
%
% NDMD( snapshots, t, Nmd )
%
% snapshots -- data matrix, each column is a time slice of data
% t -- vector of times at which data was recorded
% Nmd -- number of modes to retain

  validateattributes(t, {'numeric'},{'vector','increasing'});
  tRef = t(1);
  t = t-tRef;
  N = numel(t);
  t = t(:)'; % make into a row-vector

  if ~isscalar(Nmd)
    guess = Nmd;
    Nmd = numel(guess);
  else
    guess = complex( rand([Nmd,1]),rand([Nmd,1]) );
  end

  tM = repmat( t, [Nmd, 1] );

  dt = min(diff(t));

  validateattributes( snapshots, {'numeric'}, {'2d','ncols', numel(t)} );

  np = size(snapshots, 1); % number of observables

  [~,R] = qr( snapshots, 0 );

  function L = Lambda( z )
    z = repmat( exp(z(:)), [1,N] );
    L = z .^ tM;
    L = [L; conj(L)];
  end

  function Res = matrixRes( z )
    Res = R*( eye(N) - peye( Lambda(z) ) );
  end

  normRes = @(z) norm( matrixRes(z), 'fro').^2/2;

  disp('Guess for NDMD:')
  guess
  normRes(guess)

  lambda_a = nls_gncgs(@matrixRes, ...
                     'Jacobian-C', ...
                     complex( rand([Nmd,1]),rand([Nmd,1]) ) )
  normRes(lambda_a)

  lambda_b = ComplexOptimize( normRes, guess )
  normRes(lambda_b)

  lambda = lambda_b;

  % function v = crit( ll )
  %   ll = punpack(ll);
  %   LL = Lambda( exp(ll) );
  %   MM = R*(eye(N) - peye(LL));
  %   v = norm(MM, 'fro' );
  % end

  % l0 = complex(rand([1,Nmd]), rand([1,Nmd]));
  % opt = optimset('Diagnostics','on','Display','off','TolX',1e-12,...
  %                 'MaxFunEvals',inf, 'TolFun',1e-12,'MaxIter',2e3);

  % LB = [-inf([1,Nmd]), zeros([1,Nmd])];
  % UB = [ inf([1,Nmd]), inf([1,Nmd])];

  % lambda = punpack(fminunc( @crit, punpack(l0), opt ));
  % lambda = punpack(fmincon( @crit, punpack(l0), [],[],[],[],...
  %                           LB,UB,[],opt ));

  M = snapshots*pinv(Lambda(lambda));

end
