
close all;
clc;

%parameters
a=2;
phi=25;
l1=100;

l2=5;
% l1=5;
% l2=5;
dx=1;
dy=1;
N=10;
M=10;
ecart_type=1806;
moy=0;
%Creation of the autocorrelation function with Gaussian
x=linspace(1,N,N/dx);
y=linspace(1,M,M/dy);
[X,Y]=meshgrid(x,y);
r_gaussian=exp(-((abs((cos(phi)*X+sin(phi)*Y)/(l1)).^a)+(abs((sin(phi)*X+cos(phi)*Y)/(l2)).^a)));


%Creation of the autocorrelation function with exponential
r_exponential=exp(-((abs((cos(phi)*X+sin(phi)*Y)/(l1)).^a)+(abs((sin(phi)*X+cos(phi)*Y)/(l2)).^a)));


Spq=zeros(N,M);

for p=1:1:M
    for q=1:1:N
        n=1:N;
        m=1:M;
        Spq(p,q)=(1/M*N)*sum(sum(r_exponential(m,n).*exp(-1i*2*pi*((m*p)/M+(n*q)/N))));       
    end
end

Spq1=abs(Spq);

%creation of the normal distribute matrix and its Fourier transform
g=ecart_type*randn(N)+moy;
Gpq=zeros(N,M);
for p=1:1:M
    for q=1:1:N
        n=1:N;
        m=1:M;
        Gpq(p,q)=(1/M*N)*sum(sum(g(m,n).*exp(-1i*2*pi*((m*p)/M+(n*q)/N))));   
    end
end

figure(1);
surf(X,Y,abs(Gpq),'FaceAlpha',0.5);

%creation of the topography z
z=zeros(N,M);
for k=1:1:M
    for l=1:1:N
        s=1:N;
        r=1:M;
        z(k,l)=sum(sum(sqrt(Spq1(r,s)).*g(r,s).*exp(-1i*2*pi*((r*k)/M+(s*l)/N))));   
    end
end

full(abs(z))
figure(2);
surf(X,Y,abs(z),'FaceAlpha',0.5);