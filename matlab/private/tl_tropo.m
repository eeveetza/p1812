function Lbs = tl_tropo(dtot, theta, f, p, N0 )
%tl_tropo Basic transmission loss due to troposcatterer to P.1812-6
%   Lbs = tl_tropo(dtot, theta, f, p, N0 )
%
%   This function computes the basic transmission loss due to troposcatterer 
%   not exceeded for p% of time
%   as defined in ITU-R P.1812-6 (Section 4.4)
%
%     Input parameters:
%     dtot    -   Great-circle path distance (km)
%     theta   -   Path angular distance (mrad)
%     f       -   frequency expressed in GHz
%     p       -   percentage of time
%     N0      -   path centre sea-level surface refractivity derived 
%
%     Output parameters:
%     Lbs    -   the basic transmission loss due to troposcatterer 
%                not exceeded for p% of time
%
%     Example:
%     Lbs = tl_tropo(dtot, theta, f, p, N0 )
%       
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    01JAN16     Ivica Stevanovic, OFCOM         Initial version (P.452)
%     v1    06JUL16     Ivica Stevanovic, OFCOM         Modified version (P.1812)

%%
% Frequency dependent loss

Lf = 25*log10(f) - 2.5*(log10(f/2)).^2;    % eq (45)

% the basic transmission loss due to troposcatter not exceeded for any time
% percentage p, below 50% is given

Lbs = 190.1 + Lf + 20*log10(dtot) + 0.573*theta - 0.15*N0 - 10.125*(log10(50/p)).^(0.7);  % eq (44)

return
end
