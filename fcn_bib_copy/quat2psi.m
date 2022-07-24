function psi = quat2psi(q)

    psi = atan2( (q(2).*q(3) + q(1).*q(4) ), 0.5 - (q(3).^2 + q(4).^2) );

end