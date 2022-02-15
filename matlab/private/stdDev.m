function sigmaLoc = stdDev(f, h, R, wa)
%stdDev computes standard deviation of location variability
%   sigmaLoc = stdDev(f, h, R, wa)
%   This function computes the standard deviation according to ITU-R
%   P.1812-5, Annex 1, §4.8 and §4.10
%     Input parameters:
%     f       -   Frequency (GHz)
%     h       -   receiver/mobile antenna height above the ground (m)
%     R       -   height of representative clutter at the receiver/mobile location (m)
%     wa      -   prediction resolution, i.e., the width of the square area over which the variability applies

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    19MAR20     Ivica Stevanovic, OFCOM         Initial version
%     v1    09MAR21     Kostas Konstantinou, Ofcom      Allow h, R to be vectors

% note that there is a typo in ITU-R P.1814-5 in equation (66)
% the coefficient multiplying f should be 0.024 and not 0.0024
% this will be corrected at the 2020 WP 3K meeting

sigmaLoc = (0.52+0.024.*f) .* wa.^0.28 .* ones(size(h),class(h));

uh = ones(size(h),class(h));
IND = h<R & h>=R+10;
if any(IND)
    uh(IND) = 0;
end
IND = h<R & ~(h>=R+10);
if any(IND)
    uh(IND) = 1 - (h(IND)-R(IND))./10;
end

sigmaLoc = sigmaLoc .* uh;

return
end



