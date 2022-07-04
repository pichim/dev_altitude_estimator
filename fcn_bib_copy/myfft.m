function [X, freq] = myfft(x, Ts)

N = size(x,1);
if mod(N,2)
    N = N-1;
    x = x(1:N,:);
end

% w = hann(N, 'periodic');
% for i = 1:size(x,2)
%     x(:,i) = w.*(x(:,i) - mean(x(:,i)));
% end

freq = (0:N/2).'/(N*Ts);
X = abs(fft(x)/N);
X = X(1:N/2+1,:);
X(2:end-1,:) = 2*X(2:end-1,:);

% freq = (0:N/2).'/(N*Ts);
% X = fft(x)/N;
% P = conj(X).*X;
% P = P(1:N/2+1,:);
% P(2:end-1,:) = 2*P(2:end-1,:);
% X = sqrt(P);

end

