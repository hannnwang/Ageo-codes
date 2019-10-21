function [ ht_Cpsit ] = hwHankel(Cpsit_s,Nx,k)
%20180404
%Hankel transform applying the software package with extra format settings
%Input: Cpsit_s=Cpsit(Nx/2:end,Ny/2); k is the k for Cpsit(:,Ny/2)
%(i.e. half of k corresponds to the "Radial positions" in ht.m)
%Output: ht_Cpsit is the Hankel transform of Cpsit at all k.
k_dump=k_of_x(k);
ht_Cpsit_half=ht(Cpsit_s,k(Nx/2:end),k_dump(Nx/2:end));
ht_Cpsit=hwmakesymmetric(ht_Cpsit_half);
ht_Cpsit(Nx/2)=ht_Cpsit_half(1);
% ht_Cpsit=ht_Cpsit';

end

