function X = peye(A)
  [U,S,V] = svd(A,'econ');
  s = diag(S);
  if nargin < 2
    tol = max(size(A)) * eps(norm(s,inf));
  end
  V = V(:,s>tol);
  X = V*V';


end
