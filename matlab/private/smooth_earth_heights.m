function [hst_n, hsr_n, hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta_tot, pathtype] = smooth_earth_heights(d, h, R, htg, hrg, ae, f)
%smooth_earth_heights smooth-Earth effective antenna heights according to ITU-R P.1812-4
% [hst_n, hsr_n, hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta_tot, pathtype] = smooth_earth_heights(d, h, R, htg, hrg, ae, f)
% This function derives smooth-Earth effective antenna heights according to
% Sections 4 and 5 of the Attachment 1 to Annex 1 of ITU-R P.1812-4
%
% Input parameters:
% d         -   vector of terrain profile distances from Tx [0,dtot] (km)
% h         -   vector of terrain profile heights amsl (m)
% R         -   vector of representative clutter heights (m)
% htg, hrg  -   Tx and Rx antenna heights above ground level (m)
% ae        -   median effective Earth's radius (c.f. Eq (7a))
% f         -   frequency (GHz)
%
% Output parameters:
%
% hst_n, hsr_n -   Not corrected Tx and Rx antenna heigts of the smooth-Earth surface amsl (m)
% hst, hsr     -   Tx and Rx antenna heigts of the smooth-Earth surface amsl (m)
% hstd, hsrd   -   Tx and Rx effective antenna heigts for the diffraction model (m)
% hte, hre     -   Tx and Rx terminal effective heights for the ducting/layer reflection model (m)
% hm           -   The terrain roughness parameter (m)
% dlt          -   interfering antenna horizon distance (km)
% dlr          -   Interfered-with antenna horizon distance (km)
% theta_t      -   Interfering antenna horizon elevation angle (mrad)
% theta_r      -   Interfered-with antenna horizon elevation angle (mrad)
% theta_tot    -   Angular distance (mrad)
% pathtype     -   1 = 'los', 2 = 'transhorizon'
%
% Rev   Date        Author                          Description
% -------------------------------------------------------------------------------
% v0    15JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab (from P.452) 
% v1    15JUN16     Ivica Stevanovic, OFCOM         Modifications related to LoS path (from P.452) 
% v3    15JUN16     Ivica Stevanovic, OFCOM         Initial version for P.1812
% v4    30MAR17     Ivica Stevanovic, OFCOM         included non-corrected values of hst and hsr (87) and (88) as suggested by tranfinite
% v5    22MAR22     Ivica Stevanovic, OFCOM         updated to P.1812-6
n = length(d);

dtot = d(end);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

g = h + R;
g(1) = h(1);
g(end) = h(end);

%htc = max(hts, g(1));
%hrc = max(hrs, g(end));
htc = hts;
hrc = hrs;

% Section 5.6.1 Deriving the smooth-Earth surface

% v1 = 0;
% for ii = 2:n
%     v1 = v1 + (d(ii)-d(ii-1))*(h(ii)+h(ii-1));  % Eq (85)
% end
% v2 = 0;
% for ii = 2:n
%     v2 = v2 + (d(ii)-d(ii-1))*( h(ii)*( 2*d(ii) + d(ii-1) ) + h(ii-1) * ( d(ii) + 2*d(ii-1) ) );  % Eq (86)
% end

% the above equations optimized for speed, as suggested by Roger LeClair (leclairtelecom)

v1 = sum(diff(d) .* (h(2:n) + h(1:n-1)));  % Eq (85)
v2 = sum(diff(d) .* (h(2:n) .* (2 * d(2:n) + d(1:n-1)) + h(1:n-1) .* (d(2:n) + 2 * d(1:n-1))));  % Eq (86)

%

hst = (2*v1*dtot - v2)/dtot.^2;       % Eq (87)
hsr = (v2- v1*dtot)/dtot.^2;          % Eq (88)

hst_n = hst;
hsr_n = hsr;

% Section 5.6.2 Smooth-surface heights for the diffraction model

HH = h - (htc*(dtot-d) + hrc*d)/dtot;  % Eq (89d)

hobs = max(HH(2:n-1));                 % Eq (89a)

alpha_obt = max( HH(2:n-1)./d(2:n-1) ); % Eq (89b)

alpha_obr = max( HH(2:n-1)./( dtot - d(2:n-1) ) ); % Eq (89c)

% Calculate provisional values for the Tx and Rx smooth surface heights

gt = alpha_obt/(alpha_obt + alpha_obr);         % Eq (90e)
gr = alpha_obr/(alpha_obt + alpha_obr);         % Eq (90f)

if hobs <= 0
    hstp = hst;                                 % Eq (90a)
    hsrp = hsr;                                 % Eq (90b)
else
    hstp = hst - hobs*gt;                       % Eq (90c)
    hsrp = hsr - hobs*gr;                       % Eq (90d)
end

% calculate the final values as required by the diffraction model

if hstp >= h(1)
    hstd = h(1);                                % Eq (91a)
else
    hstd = hstp;                                % Eq (91b)
end

if hsrp > h(end)
    hsrd = h(end);                              % Eq (91c)
else
    hsrd = hsrp;                                % Eq (91d)
    
end

% Interfering antenna horizon elevation angle and distance

ii = 2:n-1;

theta =    1000 * atan( (h(ii) - hts)./(1000 * d(ii) ) - d(ii)./(2*ae) );  % Eq (77)

theta_td = 1000 * atan( (hrs   - hts)./(1000 * dtot )  -  dtot./(2*ae) );  % Eq (78)

theta_rd = 1000 * atan( (hts   - hrs)./(1000 * dtot )  -  dtot./(2*ae) );  % Eq (81)

theta_max = max(theta);                         % Eq (76)

if theta_max > theta_td                         % Eq (150): test for the trans-horizon path
    pathtype = 2; %transhorizon
else
    pathtype = 1; %los
end

theta_t = max(theta_max,theta_td);              % Eq (79)

if (pathtype == 2) %transhorizon
    
    kindex = find(theta == theta_max);
    
    lt = kindex(1)+1;   %in order to map back to path d indices, as theta takes path indices 2 to n-1, 
    
    dlt = d(lt);                                % Eq (80)
    
    % Interfered-with antenna horizon elevation angle and distance
    
    theta = 1000 * atan( (h(ii) - hrs)./(1000 * (dtot - d(ii)) ) - (dtot - d(ii))./(2*ae) );  % Eq (82a)
    
    theta_r = max(theta);
    
    kindex = find(theta == theta_r);
    
    lr = kindex(end)+1;     %in order to map back to path d indices, as theta takes path indices 2 to n-1, 
    
    dlr = dtot - d(lr);                         % Eq (83)
    
else  % pathtype == 1 (LoS)
    
    theta_r = theta_rd;                         % Eq (81)
    
    ii = 2:n-1;
    
    % speed of light as per ITU.R P.2001
    lambda = 0.2998/f;
    Ce = 1/ae;                                  % Section 4.3.1 supposing median effective Earth radius
    
    nu = (h(ii) + 500*Ce*d(ii).*(dtot-d(ii))- (hts*(dtot- d(ii)) + hrs *d(ii))/dtot).* ...
        sqrt(0.002*dtot./(lambda*d(ii).*(dtot-d(ii))));             % Eq (81)
    
    numax = max(nu);
    
    kindex = find(nu == numax);
    lt = kindex(end)+1;     %in order to map back to path d indices, as theta takes path indices 2 to n-1, 
    dlt = d(lt);                                % Eq (80)
    dlr = dtot - dlt;                           % Eq (83a)
    lr = lt;
end

% Angular distance

theta_tot = 1e3 * dtot/ae + theta_t + theta_r;  % Eq (84)


% Section 5.6.3 Ducting/layer-reflection model

% Calculate the smooth-Earth heights at transmitter and receiver as
% required for the roughness factor

hst = min(hst, h(1));                           % Eq (92a)
hsr = min(hsr, h(end));                         % Eq (92b)

% Slope of the smooth-Earth surface

m = (hsr - hst)/ dtot;                          % Eq (93)

% The terminal effective heigts for the ducting/layer-reflection model

hte = htg + h(1) -   hst;                       % Eq (94a)
hre = hrg + h(end) - hsr;                       % Eq (94b)                      

ii = lt:1:lr;

hm = max(h(ii) - (hst + m*d(ii)));              % Eq (95)

return
end