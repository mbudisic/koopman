function [DeMeaned, Mean] = removemean( Snapshots )
%REMOVEMEAN Removes time mean from snapshots.

Mean = mean(Snapshots,2);
DeMeaned = bsxfun( @minus, Snapshots, Mean );