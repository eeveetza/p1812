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
%     v3    09MAR21     Kostas Konstantinou, Ofcom      Allows R, Ct to be vectors

%%

Ah = zeros(size(R),class(R));
IND = find(h<R);
if ~isempty(IND)
    IND2 = find(Ct(IND)==3|Ct(IND)==4|Ct(IND)==5);
    if ~isempty(IND2)
        Knu = 0.342 .* sqrt(f);  % Eq (64g)
        hdif = R(IND(IND2)) - h;  % Eq (64d)
        thclut = atand(hdif./ws);  % degrees Eq (64e)
        nu = Knu .* sqrt(hdif.*thclut);  % Eq (64c)
        IND3 = nu > -0.78;
        if any(IND3)
            IND4 = IND(IND2(IND3));
            Ah(IND4) = 6.9 + 20.*log10(sqrt((nu(IND3)-0.1).^2+1)+nu(IND3)-0.1);  % Eq (12)
            Ah(IND4) = Ah(IND4) - 6.03;
        end
    end
    
    IND2 = Ct(IND)==1 | Ct(IND)==2;
    if any(IND2)
        Kh2 = 21.8 + 6.2.*log10(f);  % Eq (64f)
        Ah(IND(IND2)) = -Kh2 .* log10(h./R(IND(IND2)));
    end
    
    if any(~any(Ct(IND)==1:5,2))
        error('Unknown clutter type in function cl_loss');
    end
end
