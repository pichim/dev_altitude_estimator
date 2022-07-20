function Gds = calc_ds_frd(G, Ts, freq)

% g = squeeze(freqresp(G, 2*pi*freq));
% Warning: The frequency response is infinite at some frequencies due to poles on the stability boundary. 
g = squeeze(freqresp(G, 2*pi*freq(2:end)));
g = [g(1); g]; % copy second value to dc in case it is an integrating tf

if(isinf(g(1)))
    g(1) = g(2);
end

g(freq > 1/2/Ts) = flipud(conj(g(freq > 0 & freq < 1/2/Ts)));
Gds = frd(g, freq, Ts, 'Units', 'Hz');

end

