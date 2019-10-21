%20191021, Han Wang, hannnwangus@gmail.com
%Test the washing cycle with a test case: see if it can retrieve the corect
%W-V decomposition with the correct kappaD
%You can use this to see how figure 2 in the paper is plotted. You are also
%welcome to copy parts of this code to do your own ageostrophic anaylsis,
%but in that case, please read through the captions in this code carefully.
%I'm practising my caption writing skills, so PLEASE do not hesitate to
%contact me if you are confused about anything.

clear all
close all

%% Setting grid
kmin=7.1079e-06; %smallest k observed in data

CoarseFactor=24;%Making the grid coarser for quicker execution; setting it...
%as 1 will make k more similar to the one in data application

dk=kmin*1001/1000;
Kappas=10^(-10)*dk;
Nmax=128*24/CoarseFactor;
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
[K,L]=ndgrid(k,l);
dl=dk;

Kappa=sqrt(K.^2+L.^2);

ifig=0;
%Grids for zero padding (to be applied in the calculation of Gk)
Nmax_e=Nmax*8;
kmax=dk*Nmax_e;
k_e=linspace(-kmax+dk,kmax,2*Nmax_e);
Mmax_e=2*Nmax_e;
k_e=k_e-k_e(Mmax_e/2);
k_e=hwmakesymmetric(k_e(Mmax_e/2:end));
k_e(1:Mmax_e/2)=-k_e(1:Mmax_e/2);
Nx_e=Mmax_e;
x_e=x_of_k(k_e);
xr_e=x_e-x_e(Nx_e/2);
xr_e=hwmakesymmetric(xr_e(Mmax_e/2:end));
xr_e(1:Mmax_e/2)=-xr_e(1:Mmax_e/2);

ihalf=(Nx/2+1):Nx;

%% Some constants; they are assumed to be similar values as in data.
lat=30;%Assuming the lattitude is 30 degrees north
f= 2*7.2921*10^(-5)*sind(lat);
rho=0.36;%Assuming density is 0.36

%% Test case, where Cphit_V is calcualted from the 2D code
iks=find(k<kmin,1,'last');
Kappad_test=4*10^(-5);%$\kappa_D$ in paper
T_test=Kappad_test/2;%$T$ in paper
d=-5.3; %d and A are set just so that the magnitudes of input spectra Cut, Cvt and Cbt ...
%would somewhat look like those in data. This weired -5.3 power law does
%not mean anything and please do not think too much into it. 
A=3.07*10^(-5);
Cpsit_V_2D=A*Kappa.^d;%$S^{\psi}_V(k,l)$ in paper
Cpsit_V_2D(abs(Kappa)<kmin)=0;
Cpsit_V_2D(abs(Kappa)>(kmax-2*dk))=0;

Cphit_V_2D=...
    hwCphit_2D(Cpsit_V_2D,K,L,Lx,Ly,xr,yr,k,l,f,Kappad_test,Kappas);%$S^{\phi}_V(k,l)$ in paper
Cphit_V_2D=real(Cphit_V_2D);
Cphit_V_2D(Nx/2,Ny/2)=Cphit_V_2D(Nx/2+1,Ny/2);
Cut_V_2D=L.^2.*Cpsit_V_2D+K.^2.*Cphit_V_2D;%$S^u_V(k,l)$ in paper
Cvt_V_2D=K.^2.*Cpsit_V_2D+L.^2.*Cphit_V_2D;%$S^v_V(k,l)$ in paper
Cbt_V_2D=T_test^2*Cpsit_V_2D;%$S^b_V(k,l)$ in paper
Cut_V_k=1/(2*pi)*sum(Cut_V_2D')*dl;%$S^u_V(k)$ in paper
Cvt_V_k=1/(2*pi)*sum(Cvt_V_2D')*dl;
Cbt_V_k=1/(2*pi)*sum(Cbt_V_2D')*dl;

%wave spectra (assumed to be zero for convenience)
Cut_W_k=0.*k;%$S^u_W(k)$ in paper
Cvt_W_k=0.*k;
Cbt_W_k=0.*k;

Cut_k=Cut_W_k+Cut_V_k;
Cvt_k=Cvt_W_k+Cvt_V_k;
Cbt_k=Cbt_W_k+Cbt_V_k;

Etot_k=1/2*(Cut_k+Cvt_k+Cbt_k);%$E(k)$ in paper
EW_k_test=1/2*(Cut_W_k+Cvt_W_k+Cbt_W_k);%$E_W(k)$ in paper
EV_k_test=1/2*(Cut_V_k+Cvt_V_k+Cbt_V_k);%$E_V(k)$ in paper

%This paragraph only works when the wave spectra is set zero.
Kpsi_k_test=1/(2*pi)*sum((Kappa.^2.*Cpsit_V_2D/2)')*dl;%$K^{\psi}(k)$ in paper
Kphi_k_test=1/(2*pi)*sum((Kappa.^2.*Cphit_V_2D/2)')*dl;
EV_phi_k_test=Kphi_k_test;
EV_psi_k_test=EV_k_test-Kphi_k_test;
Kphi_k_V_test=Kphi_k_test;

ifig=ifig+1;
figure(ifig)
subplot(1,2,1)
loglog(k(ihalf)/(2*pi),EV_psi_k_test(ihalf)*rho)
hold on
loglog(k(ihalf)/(2*pi),EV_phi_k_test(ihalf)*rho)
legend('EVpsi,test case','EVphi, test case')
title('test case: 1D spectra of EVpsi and EVphi')
subplot(1,2,2)
loglog(k(ihalf)/(2*pi),Cut_k(ihalf)*rho)
hold on
loglog(k(ihalf)/(2*pi),Cvt_k(ihalf)*rho)
hold on
loglog(k(ihalf)/(2*pi),Cbt_k(ihalf)*rho)
legend('Cut','Cvt','Cbt')
xl=xlabel(['inverse wavelength $\left[m^{-1}\right]$']);
set(xl,'Interpreter','Latex','Fontsize',15);
yl=ylabel(['spectral density $\left[kg s^{-2}\right]$']);
set(yl,'Interpreter','Latex','Fontsize',15);
hold off

title('Synthetic spectra: Cut, Cvt, Cbt that ARE density-weighted')

%% Now we have finished contriving the test case. 

%One of the inputs for the decomposition in this synthetic test case is the exact
%Kphi_k. We assume in this synthetic example that we already know Kphi_k
%and Kpsi_k exactly.
Kphi_k=Kphi_k_test;
Kpsi_k=Kpsi_k_test;
%One can also reconstruct Kphi_k and Kpsi_k from Cut_k and Cvt_k applying
%the algorithm from Buhler et. al 2014. We have tried that too, and it does
%not change the results significantly in this particular test case. For
%simplicity, in this demo code, we will skip the Helmholtz decomposition
%step. Note that if we set parameters such as d or A differently and
%introduce a big difference in the magnitude between Kphi_k and Kpsi_k, say
%||Kphi_k||~0.001||Kpsi_k||,
%then the Buhler et. al. 2014 method would not correctly reconstruct Kphi_k
%because the integration error will dominate in Kphi_k. This phenomenon has
%also been noted in Person et.al 2018. 
%If you are applying this code to data, you will need to do the Helmholtz decomposition.
%The Buhler et. al. 2014 Helmholtz decomposition in this code can be
%written in this way:
% % khalf=k(ihalf);
% % Cuthalf=Cut_k(ihalf);
% % Cvthalf=Cvt_k(ihalf);
% % Cbthalf=Cbt_k(ihalf);
% % 
% % Kpsihalf=Cvthalf/2;
% % Kphihalf=Cuthalf/2;
% % for i=1:(length(khalf)-1)
% %     Kpsihalf(i)=Kpsihalf(i)+1/2/khalf(i)*trapz(khalf(i:end),Cvthalf(i:end)-Cuthalf(i:end));
% %     Kphihalf(i)=Kphihalf(i)-1/2/khalf(i)*trapz(khalf(i:end),Cvthalf(i:end)-Cuthalf(i:end));
% % end
% % 
% % %Even extension of Kpsi, Kphi
% % temp=Kpsihalf;
% % temp(end+1)=0;temp(2:end)=temp(1:end-1);
% % Kpsi_k=hwmakesymmetric(temp);
% % Kpsi_k(Nx/2)=0;
% % 
% % temp=Kphihalf;
% % temp(end+1)=0;temp(2:end)=temp(1:end-1);
% % Kphi_k=hwmakesymmetric(temp);
% % Kphi_k(Nx/2)=0;
% % 
% % %Test the Helmholtz decomp.
% % ifig=ifig+1;
% % figure(ifig)
% % subplot(1,2,1)
% % loglog(khalf,Kpsihalf)
% % hold on
% % loglog(khalf,Kpsi_k_test(ihalf),'o')
% % legend('Kpsi,BCF','Kpsi,test case')
% % subplot(1,2,2)
% % loglog(khalf,Kphihalf)
% % hold on
% % loglog(khalf,Kphi_k_test(ihalf),'o')
% % legend('Kphi,BCF','Kphi,test case')

%% Liner W-V decomposition that will be used as the initial guess for the washing cycle
%Applying the energy "equipartition" statement derived in Buhler et. al. 2014.
EW_k_BCF=2*Kphi_k;
EV_k_BCF=Etot_k-EW_k_BCF;

%% Start the washing cycle for the nonlinear decompostion
%Kappadvec(1) is smaller than the correct Kappad_test, Kappadvec(2) is
%exactly the same as Kappad_test, etc.
Kappadvec(1)=Kappad_test/2;
Kappadvec(end+1)=Kappad_test;
Kappadvec(end+1)=Kappad_test*4/3;
Kappadvec(end+1)=Kappad_test*2;

iKappa=0;
for Kappad=Kappadvec(1:end)
    iKappa=iKappa+1;
    disp('Synthetic example: washing cycle starts. The Kappa_* assumed for now is:')
    Kappad
    T=Kappad/2;%shorthand
    legendInfo{iKappa+1}=sprintf(' $\\kappa_*$=%d',double(Kappad));
    
    %initial guess for E_V^{\psi}
    EV_psi_k=EV_psi_k_test;
    
    ratio=1;
    ceps=10^(-10);
    nit=100;%maximum number of iterations
    Frelax=1;%relaxation factor (taking it to be 1 means there is no relaxation)
    iwash=0;
    ifig=ifig+1;
    fignum_dE=ifig;
    figure(fignum_dE)
    subplot(1,2,1)
    title('|dE|')
    subplot(1,2,2)
    title('KphiV')
    while(iwash<(1+nit*Frelax) && ratio>ceps &&ratio<2)
        %iwash<(1+10*Frelax) confines the control on the maximum number of
        %iterations; ratio>eps is the error control;
        %ratio>2 almost always means that the numerics is blowing up
        
        iwash=iwash+1;
        EV_psi_k_half=EV_psi_k(Nx/2+1:end);
        khalf=k(Nx/2+1:end);
        
        %Determine constant in the expression of Dpsi_V (See equation (F5) in paper)
        C0=trapz(khalf,...
            2.*khalf.*EV_psi_k_half./(T^2+khalf.^2).^(3/2));
        %Calculat Dpsi_V, which is defined as $1/(2*\pi)\int
        %l^2*S_V^{\psi}(k,l) dk$. Similar to $D^{\psi}(k)$ that
        %is introduced in Buhler et. al. 2014. 
        Dpsi_V=sqrt(khalf.^2+T^2)*C0;
        for i=2:(length(khalf))
            kco=khalf(1:i);
            ki=khalf(i);
            
            Dpsi_V(i)=Dpsi_V(i)-sqrt(ki.^2+T^2)*trapz(kco,...
                2.*kco.*EV_psi_k_half(1:i)./(T^2+kco.^2).^(3/2));
        end
        %Calculate Cpsit_V following equation (4.12) in paper
        Cpsit_V=2*EV_psi_k_half./(khalf.^2+T^2)-C0./...
            sqrt(khalf.^2+T^2);
        for i=2:(length(khalf))
            kco=khalf(1:i);
            ki=khalf(i);
            Cpsit_V(i)=Cpsit_V(i)+trapz(kco,...
                2.*kco.*EV_psi_k_half(1:i)./(T^2+kco.^2).^(3/2))./sqrt(ki.^2+T^2);
        end
        %Calculate Kpsi_V
        Kpsi_V=1/2*(khalf.^2.*Cpsit_V+Dpsi_V);
        
        %Even extension to k
        temp=Cpsit_V;
        temp(end+1)=0;temp(2:end)=temp(1:end-1);
        temp(isnan(temp))=0; %temp(isnan(temp))=temp(find(~isnan(temp), 1, 'first' ));
        p_prim=hwmakesymmetric(temp);
        p_prim(Nx/2)=p_prim(Nx/2+1);
        
        %The next steps are similar to that in the synthetic example
        %practise_Cphit_from_Cpsit_noise.m
        p=0.*k_e;
        p((Nx_e-Nx)/2+1:(Nx_e+Nx)/2)=p_prim;
        p((Nx_e-Nx)/2)=p((Nx_e+Nx)/2);
        
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
        
        G_r_spectral(Nx_e/2)=G_r_spectral(Nx_e/2+1);
        
        G_r=G_r_spectral(Nx_e/2:end);
        G_k=hwHankel(G_r',Nx_e,xr_e);
        G_k=G_k((Nx_e-Nx)/2+1:(Nx_e+Nx)/2);
        
        %testing that the subtraction trick (see Appendix E in paper, second paragraph.)
        %will not change the shape of Gk
        %significantly; we will only need this after the last iteration.
        Gkcheck=[max(G_k),mean(G_k),G_k(Nx/2),max(-(G_k-G_k(Nx/2)))];
        %the subtraction trick
        G_k=G_k-G_k(Nx/2);
        G_k(G_k<0)=0;
        
        Lhat_2D_Kphi=(Kappad^2/f./(K.^2+L.^2+Kappad^2)).^2./(K.^2+L.^2+Kappas^2)/2;
        Lhat_slice_Kphi=Lhat_2D_Kphi(:,Ny/2)';
        
        %Calculating Cphit is actually not useful in the W-V decomposition;
        %here it is calculated just for readability; readers should find it
        %easy to compare the steps to calculate Cphit here witht he steps
        %in the synthetic examples.
        Lhat_2D_Cphit=(Kappad^2/f./(K.^2+L.^2+Kappad^2)).^2./(K.^2+L.^2+Kappas^2).^2;
        Lhat_slice_Cphit=Lhat_2D_Cphit(:,Ny/2)';
        
        Kphi_slice_V=Lhat_slice_Kphi.*real(G_k);
        Kphi_slice_V(Nx/2)=Kphi_slice_V(Nx/2+1);
        Kphi_k_V=hwAbelForward(Kphi_slice_V(Nx/2:end)',Nx,k)/(2*pi);
        
        Cphit_slice=Lhat_slice_Cphit.*real(G_k);
        Cphit_slice(Nx/2)=Cphit_slice(Nx/2+1);
        Cphit_V_k=hwAbelForward(Cphit_slice(Nx/2:end)',Nx,k)/(2*pi);
        
        %New outcome for W-V decomposition
        EV_update=EV_psi_k+Kphi_k_V;
        EW_update=2*(Kphi_k-Kphi_k_V);
        
        %Update for washing cycle
        dE=Etot_k-EV_update-EW_update;
        EV_psi_k=EV_psi_k+dE/Frelax;
        
        ratio=rms(dE)/rms(Etot_k);
        
        figure(fignum_dE);
        subplot(1,2,1)
        %Checking if dE is decreasing
        loglog(k(iks+1:Nx)/(2*pi),abs(dE(iks+1:Nx))*rho);
        hold on
        subplot(1,2,2)
        loglog(k(iks+1:Nx)/(2*pi),Kphi_k_V(iks+1:Nx)*rho);
        hold on
    end
    disp('Washing cycle is terminated for this configuration. The ratio of rms between dE and E(k) is:')
    ratio
    disp('To check that the subtraction tick on G(k) is sound, these are:')
    disp('the max of G(k) before substraction, the mean of G(k) before sustraction,')
    disp('G(k) at k=0 before sustraction, and the max of the abs of negative value of G(k) after the substraction, are:')
    Gkcheck
    %have checked:
    %     Gkshow(3)/Gkshow(2) is small
    %     Gkshow(4)/Gkshow(2) is small
    disp('The real error in EV v.s. Etot is:')
    rms(EV_update-EV_k_test)/rms(Etot_k)
    
    figure(fignum_dE);
    subplot(1,2,1)
    for ilegend=1:iwash
        legend_dE{ilegend}=sprintf('|dE| at iteration %d',ilegend);
    end
    legend(legend_dE)
    hold off
    subplot(1,2,2)
    for ilegend=1:iwash
        legend_Kphi{ilegend}=sprintf('KphiV at iteration %d',ilegend);
    end
    legend_Kphi{iwash+1}=sprintf('KphiV from test case');
    loglog(k(iks+1:Nx)/(2*pi),Kphi_k_V_test(iks+1:Nx)*rho,'-.','linewidth',2.5)
    hold on
    legend(legend_Kphi,'AutoUpdate','off')
    
    %marking Kappa_d assumed in the washing cycle in the figure if Kappa_d
    %falls inside the range of k
    if Kappad>min(k(iks+1:Nx))
        xline = [Kappad/(2*pi) Kappad/(2*pi)];
        dump=ylim;
        yline = [dump(1) dump(2)];
        line(xline,yline,'Color','black','LineStyle',':');
    end
    xl=xlabel(['inverse wavelength $\left[m^{-1}\right]$']);
    set(xl,'Interpreter','Latex','Fontsize',15);
    yl=ylabel(['spectral density $\left[kg s^{-2}\right]$']);
    set(yl,'Interpreter','Latex','Fontsize',15);
    hold on
    %marking the correct Kappa_d from the test case
    xline = [Kappad_test/(2*pi) Kappad_test/(2*pi)];
    dump=ylim;
    yline = [dump(1) dump(2)];
    line(xline,yline,'Color','green','LineStyle',':');
    xl=xlabel(['inverse wavelength $\left[m^{-1}\right]$']);
    set(xl,'Interpreter','Latex','Fontsize',15);
    yl=ylabel(['spectral density $\left[kg s^{-2}\right]$']);
    set(yl,'Interpreter','Latex','Fontsize',15);
    hold off
    
    doPlotSetup
    
    %Plotting the W-V decomposition 
    ifig=ifig+1;
    figure(ifig)
    loglog(k(iks+1:Nx)/(2*pi),Etot_k(iks+1:Nx)*rho,'g')
    %The colors are not exactly the same as JAS2016 because the plots in
    %JAS2016 are probably plotted in Python
    %Do not try to specify RGB values in Matlab because if you do that,
    %Matlab cannot properly save the figures into .eps files.
    hold on
    loglog(k(iks+1:Nx)/(2*pi),EW_k_BCF(iks+1:Nx)*rho,'r');
    hold on
    loglog(k(iks+1:Nx)/(2*pi),EV_k_BCF(iks+1:Nx)*rho,'b');
    %     hold on
    %     loglog(k(iks+1:Nx)/(2*pi),EW_update(iks+1:Nx)*rho,'ro')
    hold on
    loglog(k(iks+1:Nx)/(2*pi),EV_update(iks+1:Nx)*rho,'bo')
    hold on
    loglog(k(iks+1:Nx)/(2*pi),Kphi_k_V(iks+1:Nx)*rho,'mo')
    
    %marking the correct Kappa_d from the test case
    xline = [Kappad_test/(2*pi) Kappad_test/(2*pi)];
    dump=ylim;
    yline = [dump(1) dump(2)];
    line(xline,yline,'Color','black','LineStyle',':');
    hold on
    %marking Kappa_d in the figure
    xline = [Kappad/(2*pi) Kappad/(2*pi)];
    dump=ylim;
    yline = [dump(1) dump(2)];
    line(xline,yline,'Color',[1 0.4 0],'LineStyle','-.');
    hold off
    
    leg1 = legend('$E(k)$','$E_{W}(k)$, linear','$E_{\psi,V}(k)$, linear',...
        '$E_{\psi,V}(k)$, updated',...
        '$E_{\phi,V}(k)$, updated',...
        '$\kappa_*/\left(2\pi\right)$','$\kappa_R/\left(2\pi\right)$',...
        'Location','northeast','AutoUpdate','off');
    set(leg1,'Interpreter','latex');
    
    xl=xlabel(['inverse wavelength $\left[m^{-1}\right]$']);
    set(xl,'Interpreter','Latex','Fontsize',15);
    yl=ylabel(['spectral density $\left[kg s^{-2}\right]$']);
    set(yl,'Interpreter','Latex','Fontsize',15);
    
    title('wave-vortex decomp. with ageo. correction, tropo')%this title is fashioned after JAS2016
    pbaspect([1.3 1 1])%The aspect ratio is similar to JAS2016
    
    figtitle=sprintf('tropoWV_Kappad%.2g.fig',double(Kappad));
    saveas(gcf,figtitle)
    
    savename=sprintf('tropo_WV_Kappad%.2g.mat',double(Kappad));
    save(savename)
end
