%test hwfft, hwifft


clear all
close all

% set grid points etc.
Lx=double(1e+6);%horizontal length >> water depth;
% should also be larger than the deformation radius
Ly=Lx;

Nx=512;
Ny=Nx;

x=linspace(0,Lx*(Nx-1)/Nx,Nx);
y=linspace(0,Ly*(Ny-1)/Ny,Ny);
k=k_of_x(x);
l=k_of_x(y);

xr=x-x(Nx/2);
yr=y-y(Ny/2);

%f(1:length(xr))=1;
f(1:length(xr))=cos(xr);

[ft] = hwfft(xr,k,f);
[f_v2] = hwifft(xr,k,ft);
figure(1)
subplot(3,1,1)
plot(xr,f)
subplot(3,1,2)
plot(k,ft)
subplot(3,1,3)
plot(xr,f_v2)
