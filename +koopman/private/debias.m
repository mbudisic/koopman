function [InputSnapshots, OutputSnapshots] = debias(Snapshots, r )
%DEBIAS Debias snapshots by retaining largest modes of augmented snapshots matrix.
%
%  This is the procedure from:
%    Hemati, Maziar S., and Clarence W. Rowley. 2015. “De-Biasing the Dynamic
%    Mode Decomposition for Applied Koopman Spectral Analysis.”
%    arXiv:1502.03854 [physics], February. http://arxiv.org/abs/1502.03854.
%
%  [InputSnapshots, OutputSnapshots] = debias(Snapshots)
%    Use all non-negligible SVD directions to debias Snapshots, and produce
%    Input and Output snapshot matrices.
%
%  [InputSnapshots, OutputSnapshots] = debias(Snapshots, r)
%    Retain only r dominant SVD directions to debias the matrix.
%

% We assume that OutputSnapshots = KoopmanOperator( InputSnapshots )
% column-by-column
InputSnapshots = Snapshots(:,1:end-1);
OutputSnapshots = Snapshots(:,2:end);

if nargin == 1 || ~r
  return
else
  Augmented = [InputSnapshots; OutputSnapshots];

  % Compute projector onto the non-singular right-hand space of Augmented
  % matrix
  if islogical(r) && r
    P = peye(Augmented);
  elseif isfinite(r) && r > 0
    P = peye(Augmented, ceil(r) );
  end

  InputSnapshots = InputSnapshots * P;
  OutputSnapshots = OutputSnapshots * P;
end