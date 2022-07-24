function w = getFromSkew(S)

% S = [[0 -w3 w2]; [w3 0 -w1]; [-w2 w1 0]];
w = [S(3,2), S(1,3), S(2,1)].';

end