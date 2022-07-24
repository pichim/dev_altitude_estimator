function [ang, ax, w] = quatError2AxAng(quat_des, quat)

quat_conj = quat; quat_conj(2:4) = -quat_conj(2:4);
quat_error = ( quat2QprodL(quat_conj) * quat_des.' ).';
ang = 2 * atan2( norm( quat_error(2:4) ), quat_error(1) );
ax = [1 0 0];
if abs(ang) > 1e-6
    ax = quat_error(2:4) / sin( ang / 2 );
end
w = ang * ax;

end