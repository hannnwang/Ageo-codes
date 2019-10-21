function [Pt] = hwFT2_1D_v2(Lx,xr,k,FTt1,FTt2,rflag)
%20190609
%NO ZERO PADDING to save memory
%[Pt] = hwFT2_1D(Lx,xr,k,FTt1,FTt2)
%Fourier Transform of a product of 2 terms
%PSEUDOSPECTRAL METHOD with 3/2 zero padding
%input:
%FTt1 is the FOURIER TRANSFORM of the first term
%FTt2 is the FOURIER TRANSFORM of the second term
%x,y are the horizontal and vertical grid point vectors in real space 
%Lx, Ly are horizontal and vertical domain length
%output: Pt is the Fourier Transform of the product of t1 and t2


if (nargin == 6)
    if (rflag == 1)
        % back to real space
        [A]=(hwifft(xr,k,FTt1,1));%NOTE: didn't just take the real part;
        [B]=(hwifft(xr,k,FTt2,1));
        
        % In fourier space
        Pt=hwfft(xr,k,A.*B);
    end
else
    [A]=(hwifft(xr,k,FTt1));%NOTE: didn't just take the real part;
    [B]=(hwifft(xr,k,FTt2));
    
    % In fourier space
    Pt=hwfft(xr,k,A.*B);
end

% [mem, unit] =get_free_mem()

end