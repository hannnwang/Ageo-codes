%20191021, Han Wang, hannnwangus@gmail.com
%Synthetic example: getting Cphit_V(k) from Cpsit_V(k)
%Repeat it with inputs where Cpsit_V contains noise, or where Cpsit_V is
%truncated.
%See if it works when Cpsit_V is truncated or when Cpsit_V is noisy
%Corresponding to section 3.3 in paper.
clear all
close all

dk=10^(-5);
Hs=10^(-8)*dk;
Nmax=200;
kmax=dk*Nmax;
k=linspace(-kmax+dk,kmax,2*Nmax);%length(k)=2*Nmax
Mmax=2*Nmax;%just for convenience;
k=k-k(Mmax/2);
k=hwmakesymmetric(k(Mmax/2:end));
k(1:Mmax/2)=-k(1:Mmax/2);
l=k;
Nx=Mmax;Ny=Mmax;
x=x_of_k(k);
xr=x-x(Nx/2);
xr=hwmakesymmetric(xr(Mmax/2:end));
xr(1:Mmax/2)=-xr(1:Mmax/2);
yr=xr;
Lx=2*max(xr);
Ly=Lx;
[K,L]=ndgrid(k,l);%changed to ndgrid

%Zero padding
Nmax_e=Nmax*8;
kmax=dk*Nmax_e;
k_e=linspace(-kmax+dk,kmax,2*Nmax_e);%length(k)=2*Nmax
Mmax_e=2*Nmax_e;%just for convenience;
k_e=k_e-k_e(Mmax_e/2);
k_e=hwmakesymmetric(k_e(Mmax_e/2:end));
k_e(1:Mmax_e/2)=-k_e(1:Mmax_e/2);
Nx_e=Mmax_e;
x_e=x_of_k(k_e);
xr_e=x_e-x_e(Nx_e/2);
xr_e=hwmakesymmetric(xr_e(Mmax_e/2:end));
xr_e(1:Mmax_e/2)=-xr_e(1:Mmax_e/2);

%% Test case building
%Synthetic Cpsit(k)
C=10^4;
f=7*10^(-5);
R=10000;

Cpsit=C*sqrt(pi)*R*exp(-k.^2*R^2/4);
kappaD=10^(-4);
khalf=k(Nx/2:end);

%analytical expression is known for Kphi_k_V(kappa), but the Abel forward
%has to be done numerically.
xikappa=kappaD^4*C^2.*exp(-1/8*k.^2*R^2)*pi./(f*R*(k.^2+kappaD^2)).^2;
Cphit_k_V_testcase=hwAbelForward(xikappa(Nx/2:end)',Nx,k)/(2*pi);

%% Reconstruction from "analytic" Cpsit_k
%Numerical algorithm assuming that only the projection of Cpsit is known.
p=0.*k_e;
p((Nx_e-Nx)/2+1:(Nx_e+Nx)/2)=Cpsit;
%symmetrize p
p((Nx_e-Nx)/2)=p((Nx_e+Nx)/2);
%As p is already zero at |k|>1/2|kmax|, we can skip the zero padding in
%the pseudospectral method.
G_r_spectral=1/2./(xr_e.^7).*(+32.*(xr_e).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p,(1i.*k_e).^1.*p,1),1)...
    +14.*(xr_e.^3).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^2.*p,(1i.*k_e).^2.*p,1),1)...
    -4.*(xr_e.^5).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^3.*p,(1i.*k_e).^3.*p,1),1)...
    -4.*(xr_e.^6).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^3.*p,(1i.*k_e).^4.*p,1),1)...
    -2.*(xr_e.^4).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^2.*p,(1i.*k_e).^3.*p,1),1)...
    +2.*(xr_e.^6).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^2.*p,(1i.*k_e).^5.*p,1),1)...
    -46.*(xr_e.^2).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p,(1i.*k_e).^2.*p,1),1)...
    +14.*(xr_e.^3).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p,(1i.*k_e).^3.*p,1),1)...
    -6.*(xr_e.^4).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p,(1i.*k_e).^4.*p,1),1)...
    +4.*(xr_e.^5).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p,(1i.*k_e).^5.*p,1),1)...
    +2.*(xr_e.^6).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p,(1i.*k_e).^6.*p,1),1));
%G_r blows up at 0 numerically; fix:
G_r_spectral(Nx_e/2)=G_r_spectral(Nx_e/2+1);

G_r=G_r_spectral(Nx_e/2:end);%take the x>=0 branch

G_k=hwHankel(G_r',Nx_e,xr_e);%Hankel transform
G_k=G_k((Nx_e-Nx)/2+1:(Nx_e+Nx)/2);%truncate wavenumbers

%make G_k numerically zero at k=0
G_k=G_k-G_k(Nx/2);
G_k(G_k<0)=0;%all negative values of G_k are artifacts

Lhat_2D=(kappaD^2/f./(K.^2+L.^2+kappaD^2)).^2./(K.^2+L.^2+Hs^2).^2;
Lhat_slice=Lhat_2D(:,Ny/2)';

Cphit_slice_V=Lhat_slice.*real(G_k);
Cphit_slice_V(Nx/2)=Cphit_slice_V(Nx/2+1);
Cphit_k_V_num=hwAbelForward(Cphit_slice_V(Nx/2:end)',Nx,k)/(2*pi);

%% Reconstruction from Cpsit_k "with noise"
%now, adding noise to Cpsit
Cpsit_noise = hwaddnoise_v2_1D(Cpsit,1);
p_n=0.*k_e;
p_n((Nx_e-Nx)/2+1:(Nx_e+Nx)/2)=Cpsit_noise;

p_n((Nx_e-Nx)/2)=p_n((Nx_e+Nx)/2);

G_r_spectral_n=1/2./(xr_e.^7).*(+32.*(xr_e).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_n,(1i.*k_e).^1.*p_n,1),1)...
    +14.*(xr_e.^3).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^2.*p_n,(1i.*k_e).^2.*p_n,1),1)...
    -4.*(xr_e.^5).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^3.*p_n,(1i.*k_e).^3.*p_n,1),1)...
    -4.*(xr_e.^6).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^3.*p_n,(1i.*k_e).^4.*p_n,1),1)...
    -2.*(xr_e.^4).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^2.*p_n,(1i.*k_e).^3.*p_n,1),1)...
    +2.*(xr_e.^6).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^2.*p_n,(1i.*k_e).^5.*p_n,1),1)...
    -46.*(xr_e.^2).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_n,(1i.*k_e).^2.*p_n,1),1)...
    +14.*(xr_e.^3).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_n,(1i.*k_e).^3.*p_n,1),1)...
    -6.*(xr_e.^4).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_n,(1i.*k_e).^4.*p_n,1),1)...
    +4.*(xr_e.^5).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_n,(1i.*k_e).^5.*p_n,1),1)...
    +2.*(xr_e.^6).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_n,(1i.*k_e).^6.*p_n,1),1));

G_r_spectral_n(Nx_e/2)=G_r_spectral_n(Nx_e/2+1);

G_r=G_r_spectral_n(Nx_e/2:end);

G_k_n=hwHankel(G_r',Nx_e,xr_e);
G_k_n=G_k_n((Nx_e-Nx)/2+1:(Nx_e+Nx)/2);
G_k_n=G_k_n-G_k_n(Nx/2);
G_k_n(G_k_n<0)=0;

Cphit_slice_V_n=Lhat_slice.*real(G_k_n);
Cphit_slice_V_n(Nx/2)=Cphit_slice_V_n(Nx/2+1);
Cphit_k_V_num_noise=hwAbelForward(Cphit_slice_V_n(Nx/2:end)',Nx,k)/(2*pi);

%% Reconstruction from "truncated" Cpsit_k
%now, assuming that Cpsit is not known beyond a cutoff wavenumber
kcutoff=k(Nx/2+20);
Cpsit_truncated=Cpsit;
Cpsit_truncated(abs(k)>kcutoff)=0;

p_t=0.*k_e;
p_t((Nx_e-Nx)/2+1:(Nx_e+Nx)/2)=Cpsit_truncated;

p_t((Nx_e-Nx)/2)=p_t((Nx_e+Nx)/2);

G_r_spectral_t=1/2./(xr_e.^7).*(+32.*(xr_e).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_t,(1i.*k_e).^1.*p_t,1),1)...
    +14.*(xr_e.^3).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^2.*p_t,(1i.*k_e).^2.*p_t,1),1)...
    -4.*(xr_e.^5).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^3.*p_t,(1i.*k_e).^3.*p_t,1),1)...
    -4.*(xr_e.^6).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^3.*p_t,(1i.*k_e).^4.*p_t,1),1)...
    -2.*(xr_e.^4).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^2.*p_t,(1i.*k_e).^3.*p_t,1),1)...
    +2.*(xr_e.^6).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^2.*p_t,(1i.*k_e).^5.*p_t,1),1)...
    -46.*(xr_e.^2).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_t,(1i.*k_e).^2.*p_t,1),1)...
    +14.*(xr_e.^3).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_t,(1i.*k_e).^3.*p_t,1),1)...
    -6.*(xr_e.^4).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_t,(1i.*k_e).^4.*p_t,1),1)...
    +4.*(xr_e.^5).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_t,(1i.*k_e).^5.*p_t,1),1)...
    +2.*(xr_e.^6).*hwifft(xr_e,k_e,hwFT2_1D_v2(Lx,xr_e,k_e,(1i.*k_e).^1.*p_t,(1i.*k_e).^6.*p_t,1),1));

G_r_spectral_t(Nx_e/2)=G_r_spectral_t(Nx_e/2+1);

G_r=G_r_spectral_t(Nx_e/2:end);

G_k_t=hwHankel(G_r',Nx_e,xr_e);
G_k_t=G_k_t((Nx_e-Nx)/2+1:(Nx_e+Nx)/2);
G_k_t=G_k_t-G_k_t(Nx/2);
G_k_t(G_k_t<0)=0;

Cphit_slice_V_t=Lhat_slice.*real(G_k_t);
Cphit_slice_V_t(Nx/2)=Cphit_slice_V_t(Nx/2+1);
Cphit_k_V_num_truncated=hwAbelForward(Cphit_slice_V_t(Nx/2:end)',Nx,k)/(2*pi);


fg1=figure;
iplot=(Nx/2+1):(Nx/2+round(Nx/4));%no need to plot up to the largest k 
%as the spectra decay to nearly zero there.
plot1=subplot(2,3,1);
loglog(k(iplot),Cpsit(iplot))
title('analytical')
ax1=xlabel('$k\left[m^{-1}\right]$', 'FontSize', 20);
set(ax1,'Interpreter','latex')
ax2=ylabel('$S^{\psi}_{V}(k)\left[m^5s^{-2}\right]$', 'FontSize', 20);
set(ax2,'Interpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 20)
plot2=subplot(2,3,2);
loglog(k(iplot),Cpsit_noise(iplot))
title('with noise');
ax1=xlabel('$k\left[m^{-1}\right]$', 'FontSize', 20);
set(ax1,'Interpreter','latex')
ax2=ylabel('$S^{\psi}_{V}(k)\left[m^5s^{-2}\right]$', 'FontSize', 20);
set(ax2,'Interpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 20)
plot3=subplot(2,3,3);
loglog(k(iplot),Cpsit_truncated(iplot))
title('truncated')
ax1=xlabel('$k\left[m^{-1}\right]$', 'FontSize', 20);
set(ax1,'Interpreter','latex')
ax2=ylabel('$S^{\psi}_{V}(k)\left[m^5s^{-2}\right]$', 'FontSize', 20);
set(ax2,'Interpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 20)

linkaxes([plot1,plot2,plot3],'xy')

plot4=subplot(2,3,4);
loglog(k(iplot),Cphit_k_V_num(iplot),'o')
hold on
loglog(k(iplot),Cphit_k_V_testcase(iplot))
legend('numerical','test case')
ax1=xlabel('$k\left[m^{-1}\right]$', 'FontSize', 20);
set(ax1,'Interpreter','latex')
ax2=ylabel('$S^{\phi}_{V}(k)\left[m^5s^{-2}\right]$', 'FontSize', 20);
set(ax2,'Interpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 20)
plot5=subplot(2,3,5);
loglog(k(iplot),Cphit_k_V_num_noise(iplot),'o')
hold on
loglog(k(iplot),Cphit_k_V_testcase(iplot))
legend('numerical','test case')
ax1=xlabel('$k\left[m^{-1}\right]$', 'FontSize', 20);
set(ax1,'Interpreter','latex')
ax2=ylabel('$S^{\phi}_{V}(k)\left[m^5s^{-2}\right]$', 'FontSize', 20);
set(ax2,'Interpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 20)
plot6=subplot(2,3,6);
loglog(k(iplot),Cphit_k_V_num_truncated(iplot),'o')
hold on
loglog(k(iplot),Cphit_k_V_testcase(iplot))
legend('numerical','test case')
title('truncated')
ax1=xlabel('$k\left[m^{-1}\right]$', 'FontSize', 20);
set(ax1,'Interpreter','latex')
ax2=ylabel('$S^{\phi}_{V}(k)\left[m^5s^{-2}\right]$', 'FontSize', 20);
set(ax2,'Interpreter','latex');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
yt = get(gca, 'YTick');
set(gca, 'FontSize', 20)

linkaxes([plot4,plot5,plot6],'xy')

