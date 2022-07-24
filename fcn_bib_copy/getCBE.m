function CBE = getCBE(phi, theta, psi)

% transforms earth to body or rotates a vector in body w.r.t. body frame, CBE = CEB.' = (Rz_psi*Ry_the*Rx_phi).' = Rx_phi.'*Ry_the.'*Rz_psi.'
CBE = [[                              cos(psi)*cos(theta),                              cos(theta)*sin(psi),         -sin(theta)]; ...
       [ cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(theta)*sin(phi)]; ...
       [ sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi), cos(phi)*cos(theta)]];
   
end

