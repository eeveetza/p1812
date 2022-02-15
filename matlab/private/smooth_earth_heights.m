function [hst_n, hsr_n, hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta_tot, pathtype] = smooth_earth_heights(d, h, R, htg, hrg, ae, f)
%smooth_earth_heights smooth-Earth effective antenna heights according to ITU-R P.1812-4
% [hst_n, hsr_n, hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta_tot, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f)
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
% v5    08MAR21     Kostas Konstantinou, Ofcom      d, h, R can be matrices

n = size(d,2);
dtot = d(:,end);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(:,1) + htg;
hrs = h(:,end) + hrg;
g = h + cast(R,class(h));
htc = max(hts,g(:,1));
hrc = max(hrs,g(:,end));

% Section 5.6.1 Deriving the smooth-Earth surface
tmp1 = d(:,2:n);
tmp2 = d(:,1:n-1);
tmp3 = tmp1 - tmp2;
tmp4 = h(:,1:n-1);
tmp5 = h(:,2:n);
v1 = sum(tmp3.*(tmp5+tmp4),2);  % Eq (85)
v2 = sum(tmp3.*(tmp5.*(2*tmp1+tmp2)+tmp4.*(tmp1+2*tmp2)),2);  % Eq (86)

%
hst = (2*v1.*dtot-v2) ./ dtot.^2;  % Eq (87)
hsr = (v2-v1.*dtot) ./ dtot.^2;  % Eq (88)
hst_n = hst;
hsr_n = hsr;

ii = 2 : n-1;
d2 = d(:,ii);
h2 = h(:,ii);

% Section 5.6.2 Smooth-surface heights for the diffraction model
HH = h - (htc.*(dtot-d)+hrc.*d)./dtot;  % Eq (90d)
tmp = HH(:,ii);
hobs = max(tmp,[],2);  % Eq (90a)
alpha_obt = max(tmp./d2,[],2);  % Eq (90b)
alpha_obr = max(tmp./(dtot-d2),[],2);  % Eq (90c)

% Calculate provisional values for the Tx and Rx smooth surface heights
hstp = hst;  % Eq (90a)
hsrp = hsr;  % Eq (90b)
IND = hobs > 0;
if any(IND)
    gt = alpha_obt(IND) ./ (alpha_obt(IND)+alpha_obr(IND));  % Eq (90e)
    gr = alpha_obr(IND) ./ (alpha_obt(IND)+alpha_obr(IND));  % Eq (90f)
    hstp(IND) = hst(IND) - hobs(IND).*gt;  % Eq (90c)
    hsrp(IND) = hsr(IND) - hobs(IND).*gr;  % Eq (90d)
end

% calculate the final values as required by the diffraction model
hstd = h(:,1);  % Eq (92a)
IND = hstp < h(:,1);
if any(IND)
    hstd(IND) = hstp(IND);  % Eq (91b)
end
hsrd = h(:,end);  % Eq (91c)
IND = hsrp <= h(:,end);
if any(IND)
    hsrd(IND) = hsrp(IND);  % Eq (91d)
end

% Interfering antenna horizon elevation angle and distance
theta = 1000 .* atan((h2-hts)./(1000*d2)-d2./(2*ae));  % Eq (77)
theta_td = 1000 .* atan((hrs-hts)./(1000.*dtot)-dtot./(2*ae));  % Eq (78)
theta_rd = 1000 .* atan((hts-hrs)./(1000.*dtot)-dtot./(2*ae));  % Eq (81)
theta_max = max(theta,[],2);  % Eq (76)

% Eq (150): test for the trans-horizon path
pathtype = 2 .* ones(size(theta_td),class(theta_td));  % transhorizon
IND = theta_max <= theta_td;
pathtype(IND) = 1;  % los
theta_t = max(theta_max,theta_td);              % Eq (79)
isTranshori = pathtype == 2;
lt = zeros(size(d,1),1);  % Initialise
lr = zeros(size(d,1),1);  % Initialise
dlt = zeros(size(d,1),1);  % Initialise
dlr = zeros(size(d,1),1);  % Initialise
theta_r = zeros(size(d,1),1);  % Initialise
if any(isTranshori)
    IND = isTranshori;  % Simplification
    INDf = find(IND);
    N = nnz(IND);  % Number of paths that are considered within this if statement
    dtot2 = dtot(IND);
    
    [row,col] = find(theta(IND,:)==theta_max(IND));
    if size(row,1)==1 && size(row,2)>1
        row = row.';
        col = col.';
    end
    kindex = accumarray(row,col,[N 1],@min);
    lt(IND) = kindex + 1;  % in order to map back to path d indices, as theta takes path indices 2 to n-1, 
    ind = size(d,1).*(lt(IND)-1) + INDf;
    dlt(IND) = d(ind);  % Eq (80)
    
    % Interfered-with antenna horizon elevation angle and distance
    tmp = dtot2 - d2(IND,:);
    theta = 1000 .* atan((h2(IND,:)-hrs(IND))./(1000.*tmp) - ...
        tmp./(2.*ae(IND)));  % Eq (82a)
    theta_r(IND) = max(theta,[],2);
    [row,col] = find(theta==theta_r(IND));
    if size(row,1)==1 && size(row,2)>1
        row = row.';
        col = col.';
    end
    kindex = accumarray(row,col,[N 1],@max);
    lr(IND) = kindex + 1;  % in order to map back to path d indices, as theta takes path indices 2 to n-1, 
    ind = size(d,1).*(lr(IND)-1) + INDf;
    dlr(IND) = dtot2 - d(ind);  % Eq (83)
end
if any(~isTranshori)
    IND = ~isTranshori;  % Simplification
    INDf = find(IND);
    d3 = d2(IND,:);
    dtot2 = dtot(IND);
    
    theta_r(IND) = theta_rd(IND);  % Eq (81)
    lambda = 0.3 ./ f;
    Ce = 1 ./ ae(IND);  % Section 4.3.1 supposing median effective Earth radius
    tmp = dtot2 - d3;
    nu = (h2(IND,:) + 500.*Ce.*d3.*tmp - ...
        (hts(IND).*tmp+hrs(IND).*d3)./dtot2) .* ...
        sqrt(0.002.*dtot2./(lambda.*d3.*tmp));  % Eq (81)
    numax = max(nu,[],2);
    [row,col] = find(nu==numax);
    if size(nu,1) > 1
        kindex = accumarray(row,col,[size(nu,1) 1],@max);
    else
        kindex = max(col);
    end
    lt(IND) = kindex + 1;  % in order to map back to path d indices, as theta takes path indices 2 to n-1, 
    ind = size(d,1).*(lt(IND)-1) + INDf;
    dlt(IND) = d(ind);  % Eq (80)
    dlr(IND) = dtot2 - dlt(IND);  % Eq (83a)
    lr(IND) = lt(IND);
end

% Angular distance
theta_tot = 1e3.*dtot./ae + theta_t + theta_r;  % Eq (84)

% Section 5.6.3 Ducting/layer-reflection model

% Calculate the smooth-Earth heights at transmitter and receiver as
% required for the roughness factor
hst = min(hst,h(:,1));  % Eq (92a)
hsr = min(hsr,h(:,end));  % Eq (92b)

% Slope of the smooth-Earth surface
m = (hsr-hst) ./ dtot;  % Eq (93)

% The terminal effective heigts for the ducting/layer-reflection model
hte = htg + h(:,1) - hst;  % Eq (94a)
hre = hrg + h(:,end) - hsr;  % Eq (94b)
hm = zeros(size(d,1),1,class(d));
for k = 1 : size(d,1)
    IND = lt(k) : lr(k);
    hm(k) = max(h(k,IND)-(hst(k)+m(k).*d(k,IND)));
end

return
end