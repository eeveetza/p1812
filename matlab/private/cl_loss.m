function Ah = cl_loss(h, R, Ct, f, ws)
%cl_loss additional clutter loss according to P.1812-4
%   Ah = cl_loss(h, R, Ct, f, ws)
%
%   This function computes the additional clutter loss at Tx/Rx
%   as defined in ITU-R P.1812-4 (Section 4.4)
%
%     Input parameters:
%     h       -   Antenna heigth above ground (m)
%     R       -   Clutter height (m)
%     Ct      -   Clutter type: Water/sea (1), Open/rural (2), Suburban (3)
%                 Urban/trees/forest (4), Dense urban (5)
%     f       -   frequency expressed in GHz
%     ws      -   street width parameter (default value set to 27 m)
%
%     Output parameters:
%     Ah      -   additional loss due to clutter
%
%     Example:
%      Ah = cl_loss(h, R, Ct, f, ws)
%       
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    07JUL16     Ivica Stevanovic, OFCOM         Initial version 
%     v1    17OCT16     Ivica Stevanovic, OFCOM         typo in (64f): kh2 -> Kh2 corrected
%     v2    21NOV16     Ivica Stevanovic, OFCOM         bug corrected atan -> atand

%%

if h >= R
    Ah = 0;
    return
else
    
    if Ct == 3 || Ct == 4 || Ct == 5
        
        Knu = 0.342*sqrt(f);         % Eq (64g)
        hdif = R-h;                  % Eq (64d)
        thclut = atand(hdif/ws);     % degrees Eq (64e)
        nu = Knu*sqrt(hdif*thclut);  % Eq (64c)
        
        Ah = 0;
        
        if nu > -0.78
            Ah = 6.9+ 20*log10(sqrt((nu - 0.1)^2 + 1) + nu - 0.1 );  % Eq (12)
            Ah = Ah - 6.03;
        end
        
    elseif (Ct == 1 || Ct == 2)
        Kh2 = 21.8 + 6.2*log10(f);   %Eq (64f)
        Ah = -Kh2*log10(h/R);
        
    else
        
        error('Unknown clutter type in function cl_loss');
    end
    
    return
end
        

