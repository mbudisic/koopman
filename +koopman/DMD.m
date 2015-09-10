function [lambda, Modes] = DMD(Snapshots, dt, Nmd)
  % We assume that OutputSnapshots = KoopmanOperator( InputSnapshots )
  % column-by-column
  InputSnapshots = Snapshots(:,1:end-1);
  OutputSnapshots = Snapshots(:,2:end);

  % %% De-bias input Data (Hemati, Rowley)
  % r = 20;
  % Augmented = [InputSnapshots; OutputSnapshots];
  % [Un, Sn, Vn] = svd(Augmented);
  % P = Vn(:,1:r) * Vn(:,1:r)';
  % InputSnapshots = InputSnapshots * P;
  % OutputSnapshots = OutputSnapshots * P;

  %% "Exact" DMD
  [Ux, Sx, Vx] = svd( Snapshots(:,1:end-1), 'econ' );

  % invert the Sigma matrix, skipping over almost-zero elements
  SxInv = pinv(Sx);

  %% Eigenvectors of Atilde will give Koopman modes
  Atilde = Ux' * Snapshots(:,2:end) * Vx * SxInv;

  [w, lambda] = eigs(Atilde, [], Nmd, 'LR');
  lambda = diag(lambda);

  %% Calculate modes
  %  Modes = Snapshots(:,2:end) * Vx * SxInv * w * diag(1./lambda);
  Modes = Ux * w;

  lambda = log(lambda);
  if nargin >= 2 && ~isempty(dt)
    lambda = lambda/dt;
  end
