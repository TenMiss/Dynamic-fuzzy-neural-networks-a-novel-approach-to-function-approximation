function a = RBF(d,w)
% This program is used to calculate the RBF units.
% Input:
%       d is the distance matrix
%       w is the width matrix
% Output:
%       a is the output of RBF units
% Revised 11-3-2006
% Copyright Wu Shiqian.
if nargin<2
    error('Not enough input arguments')
end
d = d .* w(:,ones(1,size(d,2)));
a = exp(-(d.*d));