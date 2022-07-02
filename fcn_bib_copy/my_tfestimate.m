function [G, C, freq, Pavg, step_resp_avg] = my_tfestimate(inp, out, window, Noverlap, N, Ts, delta)

    if(nargin == 6 || isempty(delta))
        delta = 1e-4;
    end

    % Todo: check if it is usefull to remove mean here
    inp = inp - mean(inp);
    out = out - mean(out);

    Ndata = size(inp, 1);

    % factor 2 so that the magnitude corresponds to a single sided spectrum
    % 2.3*sin(2*pi*f0*time) <=> sqrt(puu(f0)) = 2.3
    W = sum(window) / N / 2;

    Pavg = zeros(N, 3);
    Navg = 0;
    g_avg = zeros(N, 1);
    c_avg = zeros(N, 1);

    ind_start = 1;
    ind_end   = N;
    while ind_end <= Ndata

        ind = ind_start:ind_end;

        inp_act = inp(ind);
        out_act = out(ind);

        % Todo: check if it is usefull to remove mean here
        inp_act = inp_act - mean(inp_act);
        out_act = out_act - mean(out_act);

        inp_act = window .* inp_act;
        out_act = window .* out_act;

        U = fft(inp_act)/N/W;
        Y = fft(out_act)/N/W;

        Pact = [U.*conj(U) Y.*conj(U) Y.*conj(Y)];

        Pavg = Pavg + Pact;
        Navg = Navg + 1;

        [g_act, c_act] = calc_freqresp_and_cohere(Pact, delta);
        g_avg = g_avg + g_act;
        c_avg = c_avg + c_act;

        ind_start = ind_start + N - Noverlap;
        ind_end   = ind_end   + N - Noverlap;

    end

    Pavg = Pavg / Navg;
    g_avg = g_avg / Navg;
    c_avg = c_avg / Navg;
    step_resp_avg = cumsum(real(ifft(g_avg)));
    % Navg

    [g, c] = calc_freqresp_and_cohere(Pavg, delta);
    df = 1/(N*Ts);
    freq = (0:df:1/Ts-df).';

    try
        % average power spectras
        G = frd(g    , freq, Ts, 'Units', 'Hz');
        C = frd(c    , freq, Ts, 'Units', 'Hz');
    catch
        % average power spectras
        G = frd(g    , 2*pi*freq, Ts);
        C = frd(c    , 2*pi*freq, Ts);
    end
    %     % average frequency responses
    %     G = frd(g_avg, freq, Ts, 'Units', 'Hz');
    %     C = frd(c_avg, freq, Ts, 'Units', 'Hz');

end

function [g, c] = calc_freqresp_and_cohere(P, delta)

    P(:,1) = P(:,1) + delta;

    g = P(:,2) ./ P(:,1);
    c = abs(P(:,2)).^2 ./ (P(:,1) .* P(:,3));

end
