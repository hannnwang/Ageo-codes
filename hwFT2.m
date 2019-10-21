function [Pt] = hwFT2(Lx,Ly,xr,yr,k,l,FTt1,FTt2)
%Fourier Transform of a product of 2 terms
%PSEUDOSPECTRAL METHOD with 3/2 zero padding
%input:
%FTt1 is the FOURIER TRANSFORM of the first term
%FTt2 is the FOURIER TRANSFORM of the second term
%x,y are the horizontal and vertical grid point vectors in real space 
%Lx, Ly are horizontal and vertical domain length
%output: Pt is the Fourier Transform of the product of t1 and t2



[dump,Nx]=size(xr);
[dump,Ny]=size(yr);

% 3/2 zero padding method
Mx=pow2(ceil(log2(Nx+Nx/2))); My=pow2(ceil(log2(Ny+Ny/2)));

% expand x,y,k,l
xe=linspace(0,Lx*(Mx-1)/Mx,Mx);
ye=linspace(0,Ly*(My-1)/My,My);

xre=xe-xe(Mx/2);
yre=ye-ye(My/2);

ke=k_of_x(xe);
le=k_of_x(ye);

% index for the nonzero elements in the padded matrix
xstart=(Mx-Nx)/2+1;ystart=(My-Ny)/2+1;
xend=xstart+Nx-1; yend=ystart+Ny-1; 

% padding
FTt1_pad=zeros(Mx,My);
FTt1_pad(xstart:xend,ystart:yend)=FTt1;
FTt2_pad=zeros(Mx,My);
FTt2_pad(xstart:xend,ystart:yend)=FTt2;

% back to real space
[A]=(hwifft2(xre,yre,ke,le,FTt1_pad));%NOTE: didn't just take the real part; 
[B]=(hwifft2(xre,yre,ke,le,FTt2_pad));

% In fourier space
Pt_expanded=hwfft2(xre,yre,ke,le,A.*B);

% truncate the product in fourier space 
Pt=Pt_expanded(xstart:xend,ystart:yend);
end