function X = peye(A,r)
%PEYE Compute projector onto the input space of matrix A.
%
%  PEYE(A) If A is non-singular, result is the identity matrix. Otherwise, it
%  is the projector V*V' where A = USV is the economy SVD of V, where only
%  dominant singular values are retained.
%
%  PEYE(..., r) As above, except exactly r dominant singular vectors are retained.
%

% Copyright 2015 under BSD license (see LICENSE file).

  [~,S,V] = svd(A,'econ');
  s = diag(S);

  if nargin < 2
    r = nan;
  else
    try
      validateattributes(r,{'logical'},{'scalar'});
      r = nan;
    catch
      validateattributes(r,{'numeric'},...
                         {'positive','finite','integer'});
    end
  end

  if isnan(r)
    tol = max(size(A)) * eps(norm(s,inf));
    r = sum(s > tol)+1;
  end

  V(:,r:end) = [];
  s = 1./s(:);
  X = V*V';

end
