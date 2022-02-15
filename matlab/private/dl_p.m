function [Ldp,Ldb,Ld50,Lbulla50,Lbulls50,Ldsph50,maxI] = ...
    dl_p(d,g,hts,hrs,hstd,hsrd,f,omega,p,b0,DN,pol,flag4,debug)
%dl_p Diffraction loss model not exceeded for p% of time according to P.1812-4
%   function [Ldp, Ld50, Lbulla50, Lbulls50, Ldsph50] = dl_p( d, h, hts, hrs, hstd, hsrd, ap, f, omega, p, b0, DN, flag4 )
%
%   This function computes the diffraction loss not exceeded for p% of time
%   as defined in ITU-R P.1812-5 (Section 4.3-5) and Attachment 4 to Annex 1
%
%     Input parameters:
%     d       -   vector of distances di of the i-th profile point (km)
%     g       -   vector gi of heights of the i-th profile point (meters
%                 above mean sea level) + Representative clutter height. 
%                 Both vectors g and d contain n+1 profile points
%     hts     -   transmitter antenna height in meters above sea level (i=0)
%     hrs     -   receiver antenna height in meters above sea level (i=n)
%     hstd    -   Effective height of interfering antenna (m amsl) c.f. 5.6.2
%     hsrd    -   Effective height of interfered-with antenna (m amsl) c.f. 5.6.2
%     f       -   frequency expressed in GHz
%     omega   -   the fraction of the path over sea
%     p       -   percentage of time
%     b0      -   the time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100m of the lower atmosphere
%     DN      -   the average radio-refractive index lapse-rate through the
%                 lowest 1 km of the atmosphere. Note that DN is positive
%                 quantity in this procedure
%     flag4   -   Set to 1 if the alternative method is used to calculate Lbulls 
%                 without using terrain profile analysis (Attachment 4 to Annex 1)
%     pol     -   Polarization of the signal (1) horizontal, (2) vertical
 
%
%     Output parameters:
%     Ldp    -   diffraction loss for the general path not exceeded for p % of the time 
%                according to Section 4.3.4 of ITU-R P.1812-4. 
%                Ldp(1) is for the horizontal polarization 
%                Ldp(2) is for the vertical polarization
%     Ldb    -   diffraction loss for p = beta_0%
%     Ld50   -   diffraction loss for p = 50%
%     Lbulla50 -   Bullington diffraction (4.3.1) for actual terrain profile g and antenna heights
%     Lbulls50 -   Bullington diffraction (4.3.1) with all profile heights g set to zero and modified antenna heights
%     Ldshp50  -   Spherical diffraction (4.3.2) for the actual path d and modified antenna heights
%     maxI     -   Path profile index of Bullington diffraction (4.3.1) for actual terrain profile g and antenna heights
%
%
%     Example:
%     [Ldp, Ldb, Ld50, Lbulla50, Lbulls50, Ldsph50] = dl_p( d, g, hts, hrs, hstd, hsrd, ap, f, omega, p, b0, DN )
       
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    01JAN16     Ivica Stevanovic, OFCOM         Initial version (P.452)
%     v1    06JUL16     Ivica Stevanovic, OFCOM         Modifications according to P.1812
%     v2    28JUL20     Ivica Stevanovic, OFCOM         Includes Attachment 4 to Annex 1 of ITU-R P.1812-5
%                                                       with an alternative method for computation of 
%                                                       the spherical earth diffraction Lbs w/o terrain profile analysis
%     v3    09MAR21     Kostas Konstantinou, Ofcom      Allow d, g to be matrices, and hts, hrs, hstd, hsrd, omega, b0, DN to be vectors. Additional input pol.



%% 
% Use the method in 4.3.4 to calculate diffraction loss Ld for median effective 
% Earth radius ap = ae as given by equation (7a). Set median diffraction
% loss to Ldp50
[ae,ab] = earth_rad_eff(DN);
ap = ae;
[Ld50,Lbulla50,Lbulls50,Ldsph50,maxI] = ...
    dl_delta_bull(d,g,hts,hrs,hstd,hsrd,ap,f,omega,pol,flag4);

if p == 50
    Ldp = Ld50;

    if debug
        ap = ab;
        Ldb = dl_delta_bull(d,g,hts,hrs,hstd,hsrd,ap,f,omega,pol,flag4);
    else
        Ldb = [];
    end
    return
end

if p < 50
    % Use the method in 4.3.4 to calculate diffraction loss Ld for effective
    % Earth radius ap = abeta, as given in equation (7b). Set diffraction loss
    % not exceeded for beta0% time Ldb = Ld
    ap = ab;
    Ldb = dl_delta_bull(d,g,hts,hrs,hstd,hsrd,ap,f,omega,pol,flag4);

    % Compute the interpolation factor Fi
    Fi = ones(size(b0),class(b0));  % eq (40a)
    IND = p > b0;
    if any(IND)
        Fi(IND) = inv_cum_norm(p./100) ./ inv_cum_norm(b0(IND)./100);  % eq (40a)
    end
    
    % The diffraction loss Ldp not exceeded for p% of time is now given by
    Ldp = Ld50 + Fi.*(Ldb-Ld50);  % eq (41)
end

return
end