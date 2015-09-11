function norms = columnNorm( M, normfun )
%COLUMNNORM Compute norms of columns of matrix M.
%
% COLUMNNORM(M, normfun)
% Output is a row-vector, where element k is normfun(M(:,k))
% If normfun is omitted, @norm is used.

validateattributes(M, {'numeric'},{'matrix'});

if nargin == 1
    normfun = @norm;
end

norms = zeros(1,size(M,2));
for k = 1:size(M,2)
   norms(k) = normfun(M(:,k)); 
end
