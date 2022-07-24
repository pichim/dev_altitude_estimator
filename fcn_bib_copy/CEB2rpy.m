function rpy = CEB2rpy(CEB)

rpy = [atan2( CEB(3,2), CEB(3,3)), ...
       atan2(-CEB(3,1), sqrt(CEB(3,2)^2 + CEB(3,3)^2)), ...
       atan2( CEB(2,1), CEB(1,1))];
   
end

