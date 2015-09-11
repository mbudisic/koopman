function [DeMeaned, Mean] = removemean( Snapshots )
%REMOVEMEAN Removes time mean from snapshots.
%
% [Snapshots, Mean] = REMOVEMEAN(SNAPSHOTS)
%   From input data (each row in Snapshots is a time series of a single
%   observable), remove time-mean. Function returns the Snapshots file
%   without the mean, and the mean vector.

% Copyright 2015 under BSD license (see LICENSE file).

Mean = mean(Snapshots,2);
DeMeaned = bsxfun( @minus, Snapshots, Mean );