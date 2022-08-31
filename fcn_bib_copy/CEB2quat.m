function quat = CEB2quat(CEB)

% http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/

m00 = CEB(1,1); m01 = CEB(1,2); m02 = CEB(1,3);
m10 = CEB(2,1); m11 = CEB(2,2); m12 = CEB(2,3);
m20 = CEB(3,1); m21 = CEB(3,2); m22 = CEB(3,3);

tr = m00 + m11 + m22;

if (tr > 0)
  S = sqrt(tr+1.0) * 2; % S=4*qw 
  qw = 0.25 * S;
  qx = (m21 - m12) / S;
  qy = (m02 - m20) / S; 
  qz = (m10 - m01) / S; 
elseif ((m00 > m11)&(m00 > m22))
  S = sqrt(1.0 + m00 - m11 - m22) * 2; % S=4*qx 
  qw = (m21 - m12) / S;
  qx = 0.25 * S;
  qy = (m01 + m10) / S; 
  qz = (m02 + m20) / S; 
elseif (m11 > m22)
  S = sqrt(1.0 + m11 - m00 - m22) * 2; % S=4*qy
  qw = (m02 - m20) / S;
  qx = (m01 + m10) / S; 
  qy = 0.25 * S;
  qz = (m12 + m21) / S; 
else
  S = sqrt(1.0 + m22 - m00 - m11) * 2; % S=4*qz
  qw = (m10 - m01) / S;
  qx = (m02 + m20) / S;
  qy = (m12 + m21) / S;
  qz = 0.25 * S;
end

quat = [qw, qx, qy, qz];
quat = quat / norm(quat);

end

% function quat = CEB2quat(CEB)
% 
% % https://danceswithcode.net/engineeringnotes/quaternions/quaternions.html
% 
% R = CEB;
% 
% qw = sqrt( 0.25 * ( 1 + R(1,1) + R(2,2) + R(3,3) ) );
% qx = sqrt( 0.25 * ( 1 + R(1,1) - R(2,2) - R(3,3) ) );
% qy = sqrt( 0.25 * ( 1 - R(1,1) + R(2,2) - R(3,3) ) );
% qz = sqrt( 0.25 * ( 1 - R(1,1) - R(2,2) + R(3,3) ) );
% 
% quat_abs_max = max(abs([qw, qx, qy, qz]));
% 
% if quat_abs_max == abs(qw)
%     s = 0.25 / qw;
%     qx = s * (R(3,2) - R(2,3));    
%     qy = s * (R(1,3) - R(3,1));
%     qz = s * (R(2,1) - R(1,2));
% end
% 
% if quat_abs_max == abs(qx)
%     s = 0.25 / qx;
%     qw = s * (R(3,2) - R(2,3));    
%     qy = s * (R(1,2) + R(2,1));
%     qz = s * (R(1,3) + R(3,1));
% end
% 
% if quat_abs_max == abs(qy)
%     s = 0.25 / qy;
%     qw = s * (R(1,3) - R(3,1));    
%     qx = s * (R(1,2) + R(2,1));
%     qz = s * (R(2,3) + R(3,2));
% end
% 
% if quat_abs_max == abs(qz)
%     s = 0.25 / qz;
%     qw = s * (R(2,1) - R(1,2));    
%     qx = s * (R(1,3) + R(3,1));
%     qy = s * (R(2,3) + R(3,2));
% end
% 
% quat = [qw, qx, qy, qz];
% quat = quat / norm(quat);
% 
% end