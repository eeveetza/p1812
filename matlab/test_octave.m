close all
clear all
clc

if (~isOctave())
    
    fprintf(1,'This script compares results obtained in Octave\n');
    fprintf(1,'with the reference results obtained in MATLAB\n');
    fprintf(1,'using input data defined in OptimExample11Feb2022.mat\n');
    fprintf(1,'Please run the script in Octave\n');
else
    
    
    load OptimExample11Feb2022Matlab.mat
    
    [Lb_octave,Ep_octave] = tl_p1812_matr(f,p,d,h,R,Ct,zone,htg,hrg,pol,...
        'phi_t',phi_t,'phi_r',phi_r,'lam_t',lam_t,'lam_r',lam_r);
    
      
    
    delta_Lb = mean(abs(Lb - Lb_octave));
    delta_Ep = mean(abs(Ep - Ep_octave));
    
    fprintf(1,'Mean difference between Octave and Matlab Lb results: %f dB\n', delta_Lb);
    fprintf(1,'Mean difference between Octave and Matlab Ep results: %f dB\n', delta_Lb);
end

