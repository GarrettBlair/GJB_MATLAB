function b = bitslice(a,lowbit,highbit)
% BITSLICE(A,LOWBIT,HIGHBIT)
%
%   Scales return values so that the maximum is 1.0.

b = a / 2^(lowbit);
b = fix(b);
b = b / 2^(highbit - lowbit + 1);
b = b - fix(b);

b = b / max(b(:));
end