function [X, freq] = my_fft(x, Ts)

N = size(x,1);
if mod(N,2)
    N = N-1;
    x = x(1:N,:);
end

freq = (0:N/2).'/(N*Ts);
X = abs(fft(x)/N);
X = X(1:N/2+1,:);
X(2:end-1,:) = 2*X(2:end-1,:);

end

