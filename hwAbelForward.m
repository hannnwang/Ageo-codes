function [ Abel_int ] = hwAbelForward(Cpsit_s,Nx,k)
%20180522
%Switch to the numerical integral method to keep the non-negativeness

% % Previous version:
% % %20180404
% % %Abel Forward by calculating FT^{-1}H[.]
% % %The Hankel transform applys the software package with extra format settings
% % %Input: Cpsit_s=Cpsit(Nx/2:end,Ny/2); k is the k for Cpsit(:,Ny/2)
% % %(i.e. half of k corresponds to the "Radial positions" in ht.m)
% % %Output: ht_Cpsit is the Hankel transform of Cpsit at all k.
% % k_dump=k_of_x(k);
% % ht_Cpsit_half=ht(Cpsit_s,k(Nx/2:end),k_dump(Nx/2:end));
% % ht_Cpsit=hwmakesymmetric(ht_Cpsit_half);
% % ht_Cpsit(Nx/2)=ht_Cpsit_half(1);
% % Abel_ht=hwifft(k,k_dump,ht_Cpsit);
%Abel forward from the numerical integration
dk=k(2)-k(1);
Abel_int=0.*Cpsit_s;
for i=1:(length(Abel_int)-2)
    ki=k(Nx/2-1+i);
    t1=2*Cpsit_s(i)*sqrt(2*ki*dk+dk^2);
    co=k(Nx/2-1+i+1:end);
    intgd=Cpsit_s(i+1:end)'.*co./sqrt(co.^2-ki^2);
    t2=2*trapz(co,intgd);
    Abel_int(i)=t1+t2;
end
Abel_int(end-1)=2*Cpsit_s(end-1)*sqrt(2*k(end-1)*dk+dk^2);
Abel_int_sym=hwmakesymmetric(Abel_int);
Abel_int_sym(Nx/2)=Abel_int(1);
Abel_int=Abel_int_sym;
end

