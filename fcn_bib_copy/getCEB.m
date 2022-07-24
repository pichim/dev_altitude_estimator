function CEB = getCEB(phi, theta, psi)

% transforms from body to earth or rotates a vector in earth w.r.t. earth frame, CEB = Rz_psi*Ry_the*Rx_phi
CEB = [[ cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta)]; ...
       [ cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi)]; ...
       [         -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta)]];
   
end

