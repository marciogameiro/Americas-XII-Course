function [s]=cauchy_product_FFT(a1,a2)

m=length(a1);

b1=zeros(m-1,1);
ta1=[zeros(m,1);b1;a1;zeros(m,1)];
tu1=ifft(ifftshift(ta1));

b2=zeros(m-1,1);
ta2=[zeros(m,1);b2;a2;zeros(m,1)];
tu2=ifft(ifftshift(ta2));

F=fftshift(fft(tu1.*tu2));

s=real((4*m-1)*F(2*m:3*m-1));