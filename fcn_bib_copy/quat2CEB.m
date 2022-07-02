function CEB = quat2CEB(q)
% q = [qw qx qy qz], ||q|| = 1

CEB = quat2CBE(q).';

return