function [c] = time_series2Fourier(y,m)

w=ifft(y);
c=[w(end-m+2:end);w(1:m)];

end

