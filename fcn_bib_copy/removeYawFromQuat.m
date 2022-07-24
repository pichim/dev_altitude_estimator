function qn = removeYawFromQuat(q)

psi = quat2psi(q);
qz = rpy2quat([0 0 psi]);
qz(2:4) = -qz(2:4);
q = ( quat2QprodL(qz) * q.' ).';

qn = q;
n = norm(qn);
if n > 0
    qn = qn ./ n;
end

end

function psi = quat2psi(q)

    psi = atan2( (q(2).*q(3) + q(1).*q(4) ), 0.5 - (q(3).^2 + q(4).^2) );

end