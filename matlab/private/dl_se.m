function Ldsph = dl_se(d,hte,hre,ap,f,omega,pol)
%dl_se spherical-Earth diffraction loss exceeded for p% time according to ITU-R P.1812-4
%   This function computes the Spherical-Earth diffraction loss exceeded
%   for p% time for antenna heights hte and hre (m)
%   as defined in Sec. 4.3.2 of the ITU-R P.1812-4
%
%     Input parameters:
%     d       -   Great-circle path distance (km)
%     hte     -   Effective height of interfering antenna (m)
%     hre     -   Effective height of interfered-with antenna (m)
%     ap      -   the effective Earth radius in kilometers
%     f       -   frequency expressed in GHz
%     omega   -   the fraction of the path over sea
%     pol     -   Polarization of the signal (1) horizontal, (2) vertical
%
%     Output parameters:
%     Ldsph   -   The spherical-Earth diffraction loss not exceeded for p% time
%                 Ldsph(1) is for the horizontal polarization
%                 Ldsph(2) is for the vertical polarization
%
%     Example:
%     Ldsph = dl_se(d, hte, hre, ap, f, omega)

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    23DEC15     Ivica Stevanovic, OFCOM         Initial version (P.452)
%     v1    01FEB16     Ivica Stevanovic, OFCOM         Introduced dl_se_ft (P. 452)
%     v2    06JUL16     Ivica Stevanovic, OFCOM         Initial version (P.1812)
%     v3    09MAR21     Kostas Konstantinou, Ofcom      Allow d, hte, hre, ap, omega to be vectors. Additional input pol.



%% Body of function
% Wavelength in meters
lambda = 0.3/f;

% Calculate the marginal LoS distance for a smooth path
dlos = sqrt(2*ap) .* (sqrt(0.001.*hte)+sqrt(0.001.*hre));  % Eq (22)
IND = d >= dlos;
if any(IND)
    % calculate diffraction loss Ldft using the method in Sec. 4.3.3 for
    % adft = ap and set Ldsph to Ldft
    Ldsph = dl_se_ft(d,hte,hre,ap,f,omega,pol);
end
if any(~IND)
    % calculate the smallest clearance between the curved-Earth path and
    % the ray between the antennas, hse
    c = (hte-hre) ./ (hte+hre);  % Eq (24d)
    m = 250 .* d .* d ./ (ap.*(hte+hre));  % Eq (24e)
    b = 2 .* sqrt((m+1)./(3.*m)) .* ...
        cos(pi./3+1./3.*acos(3.*c./2.*sqrt(3.*m./((m+1).^3))));  % Eq (24c)
    dse1 = d ./ 2 .* (1+b);  % Eq (24a)
    dse2 = d - dse1;  % Eq (24b)
    hse = (hte-500.*dse1.*dse1./ap).*dse2 + (hre-500.*dse2.*dse2./ap).*dse1;
    hse = hse ./ d;  % Eq (23)
    
    % Calculate the required clearance for zero diffraction loss
    hreq = 17.456 .* sqrt(dse1.*dse2.*lambda./d);  % Eq (26)
    Ldsph = zeros(size(hreq),class(hreq));
    IND = find(~(hse>hreq));
    if ~isempty(IND)
        % calculate the modified effective Earth radius aem, which gives
        % marginal LoS at distance d
        aem = 500 .* (d(IND)./(sqrt(hte(IND))+sqrt(hre(IND)))).^2;  % Eq (26)
        
        % Use the method in Sec. 4.3.3 for adft = aem to obtain Ldft
        Ldft = max(dl_se_ft(d(IND),hte(IND),hre(IND),aem,f,...
            omega(IND),pol),0);
        IND2 = Ldft >= 0;
        Ldsph(IND(IND2)) = (1-hse(IND(IND2))./hreq(IND(IND2))) .* Ldft(IND2);  % Eq (27)
    end
end

return
end