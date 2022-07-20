function step_resp = calc_step_from_frd(G, C, c_min, f_max)

g = squeeze(G.ResponseData);
c = squeeze(C.ResponseData);
if isnan(abs(g(1))) % Todo: interpolate based on point 2 and 3
    g(1) = g(2);
end

g(c <= c_min) = 0;

if nargin > 3
    f = G.Frequency;
    g(f >= f_max & f <= f(end) - f_max + f(2)) = 0;
end

step_resp = cumsum(real(ifft(g)));

end