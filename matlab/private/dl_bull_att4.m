function Lbulls = dl_bull_att4(dtot, hte, hre, ap, f)
%dl_bull_att4 Bullington part of the diffraction loss for smooth path according to Attachment 4 to Annex 1 of P.1812-5
%   This function computes the spherical earth diffraction Lbuls using an
%   alternative method w/o terrain profile analysis
%   as defined in Attachment 4 to Annex 1 of ITU-R P.1812-5
%
%     Input parameters:
%     dtot    -   Great-circle path distance (km)
%     hte     -   Effective height of interfering antenna (m)
%     hre     -   Effective height of interfered-with antenna (m)
%     ap      -   the effective earth radius in kilometers
%     f       -   frequency expressed in GHz
%
%     Output parameters:
%     Lbulls   -   Bullington diffraction loss for a given path
%
%     Example:
%     Lbulls = dl_bull_att4(d, hte, hre, ap, f)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v1    28JUL20     Ivica Stevanovic, OFCOM         First implementation for P.1812-5


%% Body of function

% Effective Earth curvature Ce (km^-1)

Ce = 1/ap;

% Wavelength in meters
% speed of light as per ITU.R P.2001
lambda = 0.2998/f;


% Calculate the marginal LoS distance for a smooth path

dlos = sqrt(2*ap) * (sqrt(0.001*hte) + sqrt(0.001*hre));    % Eq (22)

if dtot < dlos % LoS
    % calculate the smallest clearance between the curved-Earth path and
    % the ray between the antennas, hse
    
    c = (hte - hre)/(hte + hre);        % Eq (24d)
    m = 250*dtot*dtot/(ap*(hte + hre));       % Eq (24e)
    
    b = 2*sqrt((m+1)/(3*m)) * cos( pi/3 + 1/3* acos( 3*c/2 * sqrt( 3*m/((m+1).^3) ) ) );   % Eq (24c)
    
    dse1 = dtot/2*(1+b);           % Eq (24a)
    dse2 = dtot - dse1;            % Eq (24b)
    
    hse = (hte - 500*dse1*dse1/ap)*dse2 + (hre - 500*dse2*dse2/ap)*dse1;
    hse = hse/dtot;                % Eq (23)
    
    % calculate the difraction parameter for the smallest clearance height hse
    % between the curved-Earth path and the ray between the antennas with the
    % distance dse1
    
    numax = -hse * sqrt(0.002*dtot/(lambda * dse1 * (dtot-dse1)) );   % Eq (105)
    
    Lus = 0;
    if numax > -0.78
        Lus = 6.9 + 20*log10(sqrt((numax-0.1).^2+1) + numax - 0.1);   % Eq (12), (106)
    end
    
else % d>=dlos, NLoS
    
    % Find the highest slope of the line from the transmitter antenna to the curved-Earth path.
    
    Stm = 500*Ce*dtot - 2*sqrt(500.0*Ce*hte);           % Eq (107)
    
    % find the highest slope of the line from the receiver antenna to the curved-Earth path
    
    Srm = 500*Ce*dtot - 2*sqrt(500.0*Ce*hre);           % Eq (108)                               % Eq (14)
    
    % Use these two slopes to calculate the Bullington point as:
    
    ds = (hre - hte + Srm*dtot)/(Stm + Srm);            % Eq (109)
    
    %Calculate the diffraction parameter nus for the Bullington point:
    
    nus = hte + Stm * ds - (hte * (dtot - ds) + hre * ds)/dtot;
    nus = nus * sqrt(0.002*dtot/(lambda * ds * (dtot - ds)));   % Eq (110)
    
    Lus = 0;
    if nus > -0.78
        Lus = 6.9 + 20*log10(sqrt((nus-0.1).^2+1) + nus - 0.1);   % Eq (12), (111)
    end
    
end
% For Luc calculated using either (106) or (111), Bullington diffraction loss
% for the path is given by

Lbulls = Lus + (1 - exp(-Lus/6.0))*(10+0.02*dtot);         % Eq (112)
return
end
