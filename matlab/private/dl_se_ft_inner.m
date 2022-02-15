function Ldft = dl_se_ft_inner(epsr,sigma,d,hte,hre,adft,f,pol)
%%dl_se_ft_inner The inner routine of the first-term spherical diffraction loss
%   This function computes the first-term part of Spherical-Earth diffraction
%   loss exceeded for p% time for antenna heights
%   as defined in Sec. 4.3.3 of the ITU-R P.1812-4, equations (29-36)
%
%     Input parameters:
%     epsr    -   Relative permittivity
%     sigma   -   Conductivity (S/m)
%     d       -   Great-circle path distance (km)
%     hte     -   Effective height of interfering antenna (m)
%     hre     -   Effective height of interfered-with antenna (m)
%     adft    -   effective Earth radius (km)
%     f       -   frequency (GHz)
%     pol     -   Polarization of the signal (1) horizontal, (2) vertical
%
%     Output parameters:
%     Ldft   -   The first-term spherical-Earth diffraction loss not exceeded for p% time
%                implementing equations (30-37), Ldft(1) is for horizontal
%                and Ldft(2) for the vertical polarization
%
%     Example:
%     Ldft = dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f)

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation (P.452)
%     v1    06JUL16     Ivica Stevanovic, OFCOM         Modifications for P.1812
%     v2    09MAR21     Kostas Konstantinou, Ofcom      Allow d, hte, hre, adft to be vectors. Additional input pol.


%% 
% For increased speed
A = (18.*sigma./f) .^ 2;

% Normalized factor for surface admittance for horizontal (1) and vertical
% (2) polarizations
K = 0.036 .* (adft.*f).^(-1/3) .* ((epsr-1).^2+A).^(-1/4);  % Eq (29a)
if pol == 2
    K = K .* (epsr.^2+A).^(1/2);  % Eq (29b)
end

% Earth ground/polarization parameter
K2 = K .^ 2;  % For increased speed
K4 = K .^ 4;  % For increased speed
beta_dft = (1+1.6.*K2+0.67.*K4) ./ (1+4.5.*K2+1.53.*K4);  % Eq (30)

% Normalized distance
X = 21.88 .* beta_dft .* (f./adft.^2).^(1/3).*d;  % Eq (31)

% Normalized transmitter and receiver heights
Yt = 0.9575 .* beta_dft .* (f.^2./adft).^(1./3).*hte;  % Eq (32a)
Yr = 0.9575 .* beta_dft .* (f.^2./adft).^(1./3).*hre;  % Eq (32b)

% Calculate the distance term given by:
IND = X >= 1.6;
Fx = zeros(size(X),class(X));
if any(IND)
    Fx(IND) = 11 + 10.*log10(X(IND)) - 17.6.*X(IND);
end
if any(~IND)
    Fx(~IND) = -20.*log10(X(~IND)) - 5.6488.*X(~IND).^1.425;  % Eq (33)
end
Bt = beta_dft .* Yt;  % Eq (35)
Br = beta_dft .* Yr;  % Eq (35)
IND = Bt > 2;
GYt = zeros(size(Bt),class(Bt));
if any(IND)
    GYt(IND) = 17.6 .* (Bt(IND)-1.1).^0.5 - 5.*log10(Bt(IND)-1.1) - 8;
end
if any(~IND)
    GYt(~IND) = 20 .* log10(Bt(~IND)+0.1.*Bt(~IND).^3);
end
IND = Br > 2;
GYr = zeros(size(Br),class(Br));
if any(IND)
    GYr(IND) = 17.6 .* (Br(IND)-1.1).^0.5 - 5.*log10(Br(IND)-1.1) - 8;
end
if any(~IND)
    GYr(~IND) = 20 .* log10(Br(~IND)+0.1.*Br(~IND).^3);
end
IND = GYr < 2+20.*log10(K);
if any(IND)
    GYr(IND) = 2 + 20.*log10(K(IND));
end
IND = GYt < 2+20.*log10(K);
if any(IND)
    GYt(IND) = 2 + 20.*log10(K(IND));
end
Ldft = -Fx - GYt - GYr;         % Eq (36)

return
end