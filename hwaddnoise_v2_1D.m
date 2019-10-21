function [Cpsit_n] = hwaddnoise_v2_1D(Cpsit,ratio)
%Add noise to Cpsit 
noise=1/ratio*(rand(size(Cpsit))-0.5);

Cpsit_n=Cpsit.*(1+noise);

end

