function [y] = buttern_filter(x0,iorder,f1,f2,dt)

% butterworth filter 
% Copyright 2016 M.Reiss and G. R?pker

% usage:
% x0 is time series, should be row
% iorder is order of butterworth filter
% f1 is lower corner frequency
% f2 is upper corner frequency
% dt is sample rate in time step

%check if x is row
if ~isrow(x0)
    x0=x0';
end

N = length(x0);
x = x0;
iswitch = 0; 

if( mod(N,2)~=0)
    x(N+1)=x0(N);
    N=N+1;
    iswitch=1;
end

% to frequncy domain
ft = fft(x);

% define cut off frequencies
omega1 = 2*pi*f1;
omega2 = 2*pi*f2;

% calculate filter coefficients
n = 2:1:N/2+1;
omega = 2*pi*(n-1)/(N*dt);
xx = (omega.*omega-omega1*omega2)./(omega.*(omega2-omega1));
prod = ones(1,length(xx));

for k=1:iorder
     prod = prod.*(xx-exp(1i*pi*(2*k-1)/(2*iorder) ) );
end

filt = (-1i)^iorder./(prod );

% apply filter
ft2 = zeros(1,N);
ft2(2:N/2+1) = ft(2:N/2+1).*filt;
ft2(N/2+2:N) = conj(fliplr(ft2(2:N/2)));
ft2 = ft2.';

% to time domain
y1 = real(ifft(ft2));

if(iswitch==1)
    y = y1(1:N-1);
else
    y = y1;
end


end
