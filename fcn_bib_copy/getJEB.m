function JEB = getJEB(phi, theta)

% jacobian of body to earth angle velocities
JEB = [[1, sin(phi)*sin(theta)/cos(theta), cos(phi)*sin(theta)/cos(theta)]; ...
       [0,                       cos(phi),                      -sin(phi)]; ...
       [0,            sin(phi)/cos(theta),            cos(phi)/cos(theta)]];
   
end

