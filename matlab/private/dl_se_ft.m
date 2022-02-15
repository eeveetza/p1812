function Ldft = dl_se_ft(d,hte,hre,adft,f,omega,pol)
%dl_se_ft First-term part of spherical-Earth diffraction according to ITU-R P.1812-4
%   This function computes the first-term part of Spherical-Earth diffraction
%   loss exceeded for p% time for antenna heights
%   as defined in Sec. 4.3.3 of the ITU-R P.1812-4
%
%     Input parameters:
%     d       -   Great-circle path distance (km)
%     hte     -   Effective height of interfering antenna (m)
%     hre     -   Effective height of interfered-with antenna (m)
%     adft    -   effective Earth radius (km)
%     f       -   frequency (GHz)
%     omega   -   fraction of the path over sea
%     pol     -   Polarization of the signal (1) horizontal, (2) vertical
%
%     Output parameters:
%     Ldft   -   The first-term spherical-Earth diffraction loss not exceeded for p% time
%                Ldft(1) is for the horizontal polarization
%                Ldft(2) is for the vertical polarization
%
%     Example:
%     Ldft = dl_se_ft(d, hte, hre, adft, f, omega)

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation (P.452)
%     v1    06JUL16     Ivica Stevanovic, OFCOM         First implementation (P.1812)
%     v2    16DEC16     Ivica Stevanovic, OFCOM         corrected bug in dl_se_ft_inner function call (epsr and sigma were interchanged)
%     v3    09MAR21     Kostas Konstantinou, Ofcom      d, hte, hre, adft, omega can be vectors. Additional input pol.


%% 
% First-term part of the spherical-Earth diffraction loss over land
IND = omega < 1;
Ldft_land = zeros(size(d),class(d));
if any(IND)
    epsr = 22;
    sigma = 0.003;
    Ldft_land(IND) = dl_se_ft_inner(epsr,sigma,d(IND),hte(IND),hre(IND),adft(IND),f,pol);
end

% First-term part of the spherical-Earth diffraction loss over sea
IND = omega > 0;
Ldft_sea = zeros(size(d),class(d));
if any(IND)
    epsr = 80;
    sigma = 5;
    Ldft_sea(IND) = dl_se_ft_inner(epsr,sigma,d(IND),hte(IND),hre(IND),adft(IND),f,pol);
end

% First-term spherical diffraction loss 
Ldft = omega.*Ldft_sea + (1-omega).*Ldft_land;  % Eq (28)

return
end
