rho=[1,0;0,0];
a=[0,0;1,0];
adagger=transpose(a);

n=80;
x = linspace(-100,100,n);
[X,P]=meshgrid(x,x);
res=eye(n);

for m=1:n;
for k=1:n;
     res(m,k) =trace(rho*expm((x(m)+1j*x(k))*adagger/2 - (x(m)-1j*x(k))*a/2 )) ; %scaled lambda
	
endfor
endfor

for m=1:2;
   for k=1,2;
      



Y=fft2(res)/(4*pi^2); %scaled lambda
%imagesc(abs(fftshift(Y)))
surf(X,P,abs(fftshift(Y)))
