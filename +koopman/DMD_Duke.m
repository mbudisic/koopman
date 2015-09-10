function [lambda, Modes] = DMD_Duke(Snapshots, dt, Nmd)
  % We assume that OutputSnapshots = KoopmanOperator( InputSnapshots )
  % column-by-column
  InputSnapshots = Snapshots(:,1:end-1);
  OutputSnapshots = Snapshots(:,2:end);

  [Q,R] = qr(InputSnapshots);

  S = pinv(R) * Q' * OutputSnapshots;

  [X, lambda] = eigs(S, [], Nmd,'LR');
  lambda = diag(lambda);

  %% Calculate modes
  Modes = InputSnapshots * X;
  lambda = log(lambda);
  if nargin >= 2 && ~isempty(dt)
    lambda = lambda/dt;
  end
