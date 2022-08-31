function Qprod = quatConj2QprodL(q)
% q = [qw qx qy qz], ||q|| = 1

q(2:4) = -q(2:4);
Qprod = quat2QprodL(q);

return