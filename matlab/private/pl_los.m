function [Lbfs, Lb0p, Lb0b] = pl_los(d,  hts, hrs, f, p, b0, dlt, dlr)
%pl_los Line-of-sight transmission loss according to ITU-R P.1812-6
%     This function computes line-of-sight transmission loss (including short-term effects)
%     as defined in ITU-R P.1812-6.
%
%     Input parameters:
%     d       -   Great-circle path distance (km)
%     hts     -   Tx antenna height above sea level (masl)
%     hrs     -   Rx antenna height above sea level (masl)
%     f       -   Frequency (GHz)
%     p       -   Required time percentage(s) for which the calculated basic
%                 transmission loss is not exceeded (%)
%     b0      -   Point incidence of anomalous propagation for the path
%                 central location (%)
%     dlt     -   For a transhorizon path, distance from the transmit antenna to
%                 its horizon (km). For a LoS path, each is set to the distance
%                 from the terminal to the profile point identified as the Bullington
%                 point in the diffraction method for 50% time
%     dlr     -   For a transhorizon path, distance from the receive antenna to
%                 its horizon (km). The same note as for dlt applies here.
%
%     Output parameters:
%     Lbfs   -   Basic transmission loss due to free-space propagation
%     Lb0p    -   Basic transmission loss not exceeded for time percentage, p%, due to LoS propagation
%     Lb0b    -   Basic transmission loss not exceedd for time percentage, b0%, due to LoS propagation
%
%     Example:
%     [Lbfs, Lb0p, Lb0b] = pl_los(d, hts, hrs, f, p, b0, dlt, dlr)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    04FEB14     Ivica Stevanovic, OFCOM         First implementation in matlab (for P.452-16)
%     v1    06JUL16     Ivica Stevanovic, OFCOM         First implementation in matlab (for P.1812-4)
%     v2    11FEB22     Ivica Stevanovic, OFCOM         Aligned with P.1812-6 to account for 3D distance and factor change


% Basic transmission loss due to free-space propagation 

%Lbfs = 92.45 + 20.0*log10(f) + 20.0*log10(d);  % (8)

dfs = sqrt(d.^2 + ((hts - hrs)/1000.0).^2);   % (8a)

Lbfs = 92.4 + 20.0*log10(f) + 20.0*log10(dfs);  % (8)

% Corrections for multipath and focusing effects at p and b0
Esp = 2.6 * (1 - exp(-0.1 * (dlt + dlr) ) ) * log10(p/50);   %(9a)
Esb = 2.6 * (1 - exp(-0.1 * (dlt + dlr) ) ) * log10(b0/50);  %(9b)

% Basic transmission loss not exceeded for time percentage p% due to
% LoS propagation
Lb0p = Lbfs + Esp;    %(10)

% Basic transmission loss not exceeded for time percentage b0% due to
% LoS propagation
Lb0b = Lbfs + Esb;    %(11)

return
end