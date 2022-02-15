function   [ae, ab] = earth_rad_eff(DN)
%earth_rad_eff Median value of the effective Earth radius
%     [ae, ab] = earth_rad_eff(DN)
%     This function computes the median value of the effective earth
%     radius, and the effective Earth radius exceeded for beta0% of time
%     as defined in ITU-R P.1812-4.
%
%     Input arguments:
%     DN      -   the average radiorefractive index lapse-rate through the
%                 lowest 1 km of the atmosphere (N-units/km)
%
%     Output arguments:
%     ae      -   the median effective Earth radius (km)
%     ab      -   the effective Earth radius exceeded for beta0 % of time
%
%     Example:
%     [ae, ab] = earth_rad_eff(DN)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    06JUL16     Ivica Stevanovic, OFCOM         First implementation in matlab
%     v1    08MAR21     Kostas Konstantinou, Ofcom      DN can be a vector


k50 = 157 ./ (157-DN);  % (6)
ae = 6371 .* k50;  % (7a)
kbeta = 3 .* ones(size(DN),class(DN));
ab = 6371 .* kbeta;  % (7b)

return
end