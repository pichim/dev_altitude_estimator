function Qprod = quat2QprodR(q)
% q = [qw qx qy qz], ||q|| = 1

% Qprod = [[ qw, -qx, -qy, -qz]; ...
%          [ qx,  qw,  qz, -qy]; ...
%          [ qy, -qz,  qw,  qx]; ...
%          [ qz,  qy, -qx,  qw]];

Qprod = [[ q(1), -q(2), -q(3), -q(4)]; ...
         [ q(2),  q(1),  q(4), -q(3)]; ...
         [ q(3), -q(4),  q(1),  q(2)]; ...
         [ q(4),  q(3), -q(2),  q(1)]];

return