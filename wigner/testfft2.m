n=16;
x = linspace(-10,10,n);
[X,P]=meshgrid(x,x);
res=cos(sqrt(X.^2 + P.^2)/2);
	
Y=fft2(res)/(4*pi^2);

%clf

imagesc(abs(fftshift(Y)))
surf(X,P,abs(fftshift(Y)))

% vacuum W
%gaus = exp(-X.^2 - P.^2)./pi;
%imagesc(gaus)
%surf(X,P,gaus)
