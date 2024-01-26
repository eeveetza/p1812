function [Lb, Ep] = tl_p1812(f, p, d, h, R, Ct, zone, htg, hrg, pol, varargin)
%tl_p1812 basic transmission loss according to P.1812-6
%   [Lb Ep] = tl_p1812(f, p, d, h, R, Ct, zone, htg, hrg, pol, varargin)
%
%   This is the MAIN function that computes the basic transmission loss not exceeded for p% time
%   and pL% locations, including additional losses due to terminal surroundings
%   and the field strength exceeded for p% time and pL% locations
%   as defined in ITU-R P.1812-6.
%   This function:
%   does not include the building entry loss (only outdoor scenarios implemented)
%
%   Other functions called from this function are in ./private/ subfolder.
%
%     Input parameters:
%     f       -   Frequency (GHz)
%     p       -   Required time percentage for which the calculated basic
%                 transmission loss is not exceeded
%     d       -   vector of distances di of the i-th profile point (km)
%     h       -   vector of heights hi of the i-th profile point (meters
%                 above mean sea level.
%     R       -   vector of representative clutter height Ri of the i-th profile point (m)
%     Ct      -   vector of representative clutter type Cti of the i-th profile point
%                 Water/sea (1), Open/rural (2), Suburban (3),
%                 Urban/trees/forest (4), Dense urban (5)
%                 if empty or all zeros, the default clutter used is Open/rural
%     zone    -   vector of radio-climatic zone types: Inland (4), Coastal land (3), or Sea (1)
%     htg     -   Tx Antenna center heigth above ground level (m)
%     hrg     -   Rx Antenna center heigth above ground level (m)
%     pol     -   polarization of the signal (1) horizontal, (2) vertical
%
%    Input parameters related to path centre:
%    EITHER the following are required:
%
%     phi_t    - latitude of Tx station (degrees)
%     phi_r    - latitude of Rx station (degrees)
%     lam_t    - longitude of Tx station (degrees)
%     lam_r    - longitude of Rx station (degrees)
%
%    OR the following are required:
%
%     phi_path - latitude of the path centre (degrees)
%
%   Examples of both cases are provided below.
%
%   Output parameters:
%     Lb   - basic transmission loss according to P.1812-6
%     Ep   - the field strength relative to Ptx
%
%
% Examples:
%
% 1) Call with required input parameters, and latitude/longitude of Tx/Rx:
% [Lb,Ep] = tl_p1812(f,p,d,h,R,Ct,zone,htg,hrg,pol,...
%     'phi_t',phi_t,'phi_r',phi_r,'lam_t',lam_t,'lam_r',lam_r)
%
% 2) Call with required input parameters, and latitude of path centre:
% [Lb,Ep] = tl_p1812(f,p,d,h,R,Ct,zone,htg,hrg,pol,'phi_path',phi_path);
%
% 3) Call with Name-Value Pair Arguments. Name is the argument name and Value is the
% corresponding value. Name must appear inside quotes.
% [Lb,Ep] = tl_p1812(___,Name,Value)
% Example: tl_p1812(f,p,d,h,R,Ct,zone,htg,hrg,pol,'phi_path',phi_path,'DN',DN,'N0',N0)
% Below are the valid Name-Value Pair Arguments:
%     pL      -   Required time percentage for which the calculated basic
%                 transmission loss is not exceeded (1% - 99%)
%     sigmaL  -   location variability standard deviations computed using
%                 stdDev.m according to §4.8 and §4.10
%                 the value of 5.5 dB used for planning Broadcasting DTT
%     Ptx     -   Transmitter power (kW), default value 1 kW
%     DN      -   The average radio-refractive index lapse-rate through the
%                 lowest 1 km of the atmosphere (it is a positive quantity in this
%                 procedure) (N-units/km)
%     N0      -   The sea-level surface refractivity, is used only by the
%                 troposcatter model as a measure of location variability of the
%                 troposcatter mechanism. The correct values of DN and N0 are given by
%                 the path-centre values as derived from the appropriate
%                 maps (N-units)
%     dct     -   Distance over land from the transmit and receive
%     dcr         antennas to the coast along the great-circle interference path (km).
%                 default values dct = 500 km, dcr = 500 km, or
%                 set to zero for a terminal on a ship or sea platform
%     flag4   -   Set to 1 if the alternative method is used to calculate Lbulls
%                 without using terrain profile analysis (Attachment 4 to Annex 1)
%     debug   -   Set to 1 if the log files are to be written,
%                 otherwise set to default 0
%     fid_log  -   if debug == 1, a file identifier of the log file can be
%                 provided, if not, the default file with a file
%                 containing a timestamp will be created
%



%
% Numbers refer to Rec. ITU-R P.1812-6

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    02JAN16     Ivica Stevanovic, OFCOM         Initial version for ITU-R P.452-16
%     v1    07JUL16     Ivica Stevanovic, OFCOM         Modified for ITU-R P.1812-4
%     v2    23AUG16     Ivica Stevanovic, OFCOM         Included path
%                                                       center calculation
%     v3    17OCT16     Ivica Stevanovic, OFCOM         typo in cl_loss call corrected (J.Dieterle)
%                                                       h(1)-> htg, h(end)-> hrg
%     v4    03MAY18     Ivica Stevanovic, OFCOM         Introduced path center calculation
%                                                       according to Annex
%                                                       H of ITU-R P.2001-2
%                                                       and additional validation checks
%     v5    28MAY18     Ivica Stevanovic, OFCOM         Corrected a bug in printing out Lbulla50 and Lbulls50
%                                                       (scalar and not vector values), reported by Damian Bevan
%     v6    08JUN18     Ivica Stevanovic, OFCOM         Corrected a bug (swap of long and lat coordinates in
%                                                       great_circle_path call, reported by Wladislaw Budynkiewicz.
%     v7    19MAR20     Ivica Stevanovic, OFCOM         Aligned with P.1812-5, introduced location variability
%     v8    28JUL20     Ivica Stevanovic, OFCOM         Introduced alternative method to compute Lbulls w/o using
%                                                       terrain profile (Attachment 4 to Annex 1)
%     v9    08MAR21     Kostas Konstantinou, Ofcom      Simplified handling of input arguments using name-value pairs
%     v10   11FEB22     Ivica Stevanovic, OFCOM         Aligned with P.1812-6, renamed subfolder "src" into "private"
%                                                       which is automatically in MATLAB search path ..
%     v11   10MAR22     Ivica Stevanovic, OFCOM         Use comma as a separator in the written csv files instead of semicolon
%

% MATLAB Version 9.2.0.556344 (R2017a) used in development of this code
%
% The Software is provided "AS IS" WITH NO WARRANTIES, EXPRESS OR IMPLIED, 
% INCLUDING BUT NOT LIMITED TO, THE WARRANTIES OF MERCHANTABILITY, FITNESS 
% FOR A PARTICULAR PURPOSE AND NON-INFRINGMENT OF INTELLECTUAL PROPERTY RIGHTS 
% 
% Neither the Software Copyright Holder (or its affiliates) nor the ITU 
% shall be held liable in any event for any damages whatsoever
% (including, without limitation, damages for loss of profits, business 
% interruption, loss of information, or any other pecuniary loss)
% arising out of or related to the use of or inability to use the Software.
%
% THE AUTHOR(S) AND OFCOM (CH) DO NOT PROVIDE ANY SUPPORT FOR THIS SOFTWARE
%
% This function calls other functions that are placed in the ./private folder


%% Read the input arguments and check them

% Parse the optional input
iP = inputParser;
iP.addParameter('phi_path',[],@(x) isnumeric(x) && abs(x)<=80 )
iP.addParameter('phi_t',[],@(x) isnumeric(x) && abs(x)<=80)
iP.addParameter('phi_r',[],@(x) isnumeric(x) && abs(x)<=80)
iP.addParameter('lam_t',[],@(x) isnumeric(x))
iP.addParameter('lam_r',[],@(x) isnumeric(x))
iP.addParameter('pL',50,@(x) isnumeric(x) && x>=1 && x<=99)
iP.addParameter('sigmaL',0)
iP.addParameter('Ptx',1)
iP.addParameter('DN',45)
iP.addParameter('N0',325)
iP.addParameter('dct',500)
iP.addParameter('dcr',500)
iP.addParameter('flag4',0)
iP.addParameter('debug',0)
iP.addParameter('fid_log',[])
iP.parse(varargin{:});

% Unpack from input parser
dcr = iP.Results.dcr;
dct = iP.Results.dct;
debug = iP.Results.debug;
DN = iP.Results.DN;
fid_log = iP.Results.fid_log;
flag4 = iP.Results.flag4;
N0 = iP.Results.N0;
phi_path = iP.Results.phi_path;
if isempty(phi_path)
    phi_t = iP.Results.phi_t;
    phi_r = iP.Results.phi_r;
    lam_t = iP.Results.lam_t;
    lam_r = iP.Results.lam_r;
    if(~isOctave())
        mustBeNonempty(phi_t)
        mustBeNonempty(phi_r)
        mustBeNonempty(lam_t)
        mustBeNonempty(lam_r)
    end
end
pL = iP.Results.pL;
Ptx = iP.Results.Ptx;
sigmaL = iP.Results.sigmaL;


% verify input argument values and limits
check_limit(f, 0.03, 6.0, 'f [GHz]');
check_limit(p, 1, 50, 'p [%]');
check_limit(htg, 1, 3000, 'htg [m]');
check_limit(hrg, 1, 3000, 'hrg [m]');
check_value(pol, [1, 2], 'Polarization (pol) ');
%check_value(Ct, [1, 2, 3, 4, 5], 'Clutter coverage (Ct) ');
check_value(zone, [1, 3, 4], 'Radio-climatic zone (zone) ');

NN=length(d);

% the number of elements in d and path need to be the same
if(length(h) ~= NN)
    error('The number of elements in the array ''d'' and array ''h'' must be the same.')
end

if isempty(R)
    R=zeros(size(h)); % default is clutter height zero
else
    if(length(R) ~= NN)
        error('The number of elements in the array ''d'' and array ''R'' must be the same.')
    end
end

if isempty(Ct)
    Ct = 2*ones(size(h)); % default is Open/rural clutter type
    
elseif Ct == 0
    Ct = 2*ones(size(h)); % default is Open/rural clutter type
else
    if(length(Ct) ~= NN)
        error('The number of elements in the array ''d'' and array ''Ct'' must be the same.')
    end
end

if isempty(zone)
    zone = 4*ones(size(h)); % default is Inland radio-meteorological zone
else
    if(length(zone) ~= NN)
        error('The number of elements in the array ''d'' and array ''zone'' must be the same.')
    end
end

if zone(1) == 1 % Tx at sea
    dct = 0;
end
if zone(end) == 1 %Rx at sea
    dcr = 0;
end


% handle number fidlog is reserved here for writing the files
% if fidlog is already open outside of this function, the file needs to be
% empty (nothing written), otherwise it will be closed and opened again
% if fid is not open, then a file with a name corresponding to the time
% stamp will be opened and closed within this function.

inside_file = 0;
if (debug)
    if (isempty (fid_log))
        fid_log = fopen(['P1812_' num2str(floor(now*1e10)) '_log.csv'],'w');
        inside_file = 1;
        if (fid_log == -1)
            error('The log file could not be opened.')
        end
    else
        inside_file = 0;
    end
end

floatformat= '%.10g,\n';

if (debug)
    
    fprintf(fid_log,'# Parameter,Ref,,Value,\n');
    fprintf(fid_log,['Ptx (kW),,,' floatformat],Ptx);
    fprintf(fid_log,['f (GHz),,,' floatformat],f);
    fprintf(fid_log,['p (%%),,,' floatformat],p);
    fprintf(fid_log,['pL (%%),,,' floatformat],pL);
    fprintf(fid_log,['sigmaL (dB),,,' floatformat],sigmaL);
    fprintf(fid_log,['phi_t (deg),,,' floatformat],phi_t);
    fprintf(fid_log,['phi_r (deg),,,' floatformat],phi_r);
    fprintf(fid_log,['lam_t (deg),,,' floatformat],lam_t);
    fprintf(fid_log,['lam_r (deg),,,' floatformat],lam_r);
    fprintf(fid_log,['htg (m),,,' floatformat],htg);
    fprintf(fid_log,['hrg (m),,,' floatformat],hrg);
    fprintf(fid_log,['pol,,,' '%d,\n'],pol);
    %    fprintf(fid_log,['ws (m);;;' floatformat],ws);
    fprintf(fid_log,['DN ,,,' floatformat],DN);
    fprintf(fid_log,['N0 ,,,' floatformat],N0);
    fprintf(fid_log,['dct (km) ,,,' floatformat],dct);
    fprintf(fid_log,['dcr (km) ,,,' floatformat],dcr);
    fprintf(fid_log,['R2 (m) ,,,' floatformat],R(2));
    fprintf(fid_log,['Rn-1 (m) ,,,' floatformat],R(end-1));
    fprintf(fid_log,['Ct Tx ,Table 2,,' floatformat],Ct(2));
    fprintf(fid_log,['Ct Rx ,Table 2,,' floatformat],Ct(end-1));
    
end


% Compute the path profile parameters

% Path center latitude
if isempty(phi_path)
    Re = 6371;
    dpnt = 0.5 .* (d(end)-d(1));
    [~,phi_path] = great_circle_path(lam_r,lam_t,phi_r,phi_t,Re,dpnt);
end

% Compute  dtm     -   the longest continuous land (inland + coastal =34) section of the great-circle path (km)
zone_r = 34;
dtm = longest_cont_dist(d, zone, zone_r);

% Compute  dlm     -   the longest continuous inland section (4) of the great-circle path (km)
zone_r = 4;
dlm = longest_cont_dist(d, zone, zone_r);

% Compute b0
b0 = beta0(phi_path, dtm, dlm);

[ae, ab] = earth_rad_eff(DN);

% Compute the path fraction over sea Eq (1)

omega = path_fraction(d, zone, 1);

% Derive parameters for the path profile analysis

[hst_n, hsr_n, hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, R, htg, hrg, ae, f);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Modify the path by adding representative clutter, according to Section 3.2
% excluding the first and the last point
g = h + R;
g(1) = h(1);
g(end) = h(end);

%Compute htc and hrc as defined in Table 5 (P.1812-6)
% htc = max(hts,g(1));
% hrc = max(hrs,g(end));
htc = hts;
hrc = hrs;

if (debug)
    fprintf(fid_log,',,,,\n');
    %fprintf(fid_log,'# Path;Ref;;Value;\n');
    fprintf(fid_log,['d (km),,,' floatformat],dtot);
    fprintf(fid_log,['dlt (km),Eq (78),,' floatformat],dlt);
    fprintf(fid_log,['dlr (km),Eq (81a),,' floatformat],dlr);
    fprintf(fid_log,['th_t (mrad),Eqs (76-78),,' floatformat],theta_t);
    fprintf(fid_log,['th_r (mrad),Eqs (79-81),,' floatformat],theta_r);
    fprintf(fid_log,['th (mrad),Eq (82),,' floatformat],theta);
    fprintf(fid_log,['hts (m),,,' floatformat],hts);
    fprintf(fid_log,['hrs (m),,,' floatformat],hrs);
    fprintf(fid_log,['htc (m),Table 5,,' floatformat],htc);
    fprintf(fid_log,['hrc (m),Table 5,,' floatformat],hrc);
    fprintf(fid_log,['w,Table 5,,' floatformat],omega);
    fprintf(fid_log,['dtm (km),Sec 3.6,,' floatformat],dtm);
    fprintf(fid_log,['dlm (km),Sec 3.6,,' floatformat],dlm);
    fprintf(fid_log,['phi (deg),Eq (4),,' floatformat],phi_path);
    fprintf(fid_log,['b0 (%%),Eq (5),,' floatformat],b0);
    fprintf(fid_log,['ae (km),Eq (7a),,' floatformat],ae);
    fprintf(fid_log,['hst (m),Eq (85),,' floatformat],hst_n);
    fprintf(fid_log,['hsr (m),Eq (86),,' floatformat],hsr_n);
    fprintf(fid_log,['hst (m),Eq (90a),,' floatformat],hst);
    fprintf(fid_log,['hsr (m),Eq (90b),,' floatformat],hsr);
    fprintf(fid_log,['hstd (m),Eq (89),,' floatformat],hstd);
    fprintf(fid_log,['hsrd (m),Eq (89),,' floatformat],hsrd);
    fprintf(fid_log,['htc'' (m),Eq (37a),,' floatformat],htc-hstd);
    fprintf(fid_log,['hrc'' (m),Eq (37b),,' floatformat],hrc-hsrd);
    fprintf(fid_log,['hte (m),Eq (92a),,' floatformat],hte);
    fprintf(fid_log,['hre (m),Eq (92b),,' floatformat],hre);
    fprintf(fid_log,['hm (m),Eq (93),,' floatformat],hm);
    fprintf(fid_log,'\n');
    
    
end


% Calculate an interpolation factor Fj to take account of the path angular
% distance Eq (57)

THETA = 0.3;
KSI = 0.8;

Fj = 1.0 - 0.5*( 1.0 + tanh(3.0 * KSI * (theta-THETA)/THETA) );

% Calculate an interpolation factor, Fk, to take account of the great
% circle path distance:

dsw = 20;
kappa = 0.5;

Fk = 1.0 - 0.5*( 1.0 + tanh(3.0 * kappa * (dtot-dsw)/dsw) ); % eq (58)

%[Lbfs, Lb0p, Lb0b] = pl_los(dtot, f, p, b0, dlt, dlr);
[Lbfs, Lb0p, Lb0b] = pl_los(dtot, hts, hrs, f, p, b0, dlt, dlr);

[ Ldp, Ldb, Ld50, Lbulla50, Lbulls50, Ldsph50] = dl_p( d, g, htc, hrc, hstd, hsrd, f, omega, p, b0, DN, flag4 );

% The median basic transmission loss associated with diffraction Eq (42)

Lbd50 = Lbfs + Ld50;

% The basic tranmission loss associated with diffraction not exceeded for
% p% time Eq (43)

Lbd = Lb0p + Ldp;

% A notional minimum basic transmission loss associated with LoS
% propagation and over-sea sub-path diffraction

Lminb0p = Lb0p + (1-omega)*Ldp;

% eq (40a)
Fi = 1;

if p >= b0
    
    Fi = inv_cum_norm(p/100)/inv_cum_norm(b0/100);
    
    Lminb0p = Lbd50 + (Lb0b + (1-omega)*Ldp - Lbd50)*Fi;   % eq (59)
    
end


% Calculate a notional minimum basic transmission loss associated with LoS
% and transhorizon signal enhancements

eta = 2.5;

Lba = tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, omega, ae, b0);

Lminbap = eta*log(exp(Lba/eta) + exp(Lb0p/eta));           % eq (60)

% Calculate a notional basic transmission loss associated with diffraction
% and LoS or ducting/layer reflection enhancements

Lbda = Lbd;

if Lminbap <= Lbd(1)
    Lbda(1) = Lminbap + (Lbd(1)-Lminbap)*Fk; % eq (61)
end

if Lminbap <= Lbd(2)
    Lbda(2) = Lminbap + (Lbd(2)-Lminbap)*Fk; % eq (61)
end

% Calculate a modified basic transmission loss, which takes diffraction and
% LoS or ducting/layer-reflection enhancements into account

Lbam = Lbda + (Lminb0p - Lbda)*Fj;   % eq (62)

% Calculate the basic transmission loss due to troposcatter not exceeded
% for any time percantage p

Lbs = tl_tropo(dtot, theta, f, p, N0);

% Calculate the final transmission loss not exceeded for p% time
% ignoring the effects of terminal clutter

Lbc_pol = -5*log10(10.^(-0.2*Lbs) + 10.^(-0.2*Lbam)) ;  % eq (63)

% choose the right polarization

Lbc = Lbc_pol(pol);

% % The additional clutter losses from are removed in P.1812-6
% % additional losses due to terminal surroundings (Section 4.7)
%
% % Parameter ws relates to the width of the street. It is set to 27 unless
% % there is specific local information available
%
% ws = 27;
%
% % Transmitter side
%
% Aht = cl_loss(htg, R(1), Ct(1), f, ws);
%
% % Receiver side
%
% Ahr = cl_loss(hrg, R(end), Ct(end), f, ws);
%
% % Basic transmission loss not exceeded for p% time and 50% locations,
% % including the effects of terminal clutter losses

%Lbc = Lbu + Aht + Ahr;

Lloc = 0.0; % outdoors only (67a)

% Location variability of losses (Section 4.8)
if zone(end) ~= 1 %Rx at sea
    Lloc = -inv_cum_norm( pL/100 ) * sigmaL;
end

% Basic transmission loss not exceeded for p% time and pL% locations
% (Sections 4.8 and 4.9) not implemented

Lb = max(Lb0p, Lbc + Lloc);  %  eq (69)

% The field strength exceeded for p% time and pL% locations

Ep = 199.36+ 20*log10(f) - Lb; % eq (70)

% % Scale to the transmitter power

EpPtx = Ep + 10*log10(Ptx);

if (debug)
    fprintf(fid_log,['Fi,Eq (40),,' floatformat],Fi);
    fprintf(fid_log,['Fj,Eq (57),,' floatformat],Fj);
    fprintf(fid_log,['Fk,Eq (58),,' floatformat],Fk);
    fprintf(fid_log,['Lbfs,Eq (8),,' floatformat],Lbfs);
    fprintf(fid_log,['Lb0p,Eq (10),,' floatformat],Lb0p);
    fprintf(fid_log,['Lb0b,Eq (11),,' floatformat],Lb0b);
    fprintf(fid_log,['Lbulla (dB),Eq (21),,' floatformat],Lbulla50);
    fprintf(fid_log,['Lbulls (dB),Eq (21),,' floatformat],Lbulls50);
    fprintf(fid_log,['Ldsph (dB),Eq (27),,' floatformat],Ldsph50(pol));
    fprintf(fid_log,['Ld50 (dB),Eq (39),,' floatformat],Ld50(pol));
    fprintf(fid_log,['Ldb (dB),Eq (39),,' floatformat],Ldb(pol));
    fprintf(fid_log,['Ldp (dB),Eq (41),,' floatformat],Ldp(pol));
    fprintf(fid_log,['Lbd50 (dB),Eq (42),,' floatformat],Lbd50(pol));
    fprintf(fid_log,['Lbd (dB),Eq (43),,' floatformat],Lbd(pol));
    
    fprintf(fid_log,['Lminb0p (dB),Eq (59),,' floatformat],Lminb0p(pol));
    
    fprintf(fid_log,['Lba (dB),Eq (46),,' floatformat],Lba);
    fprintf(fid_log,['Lminbap (dB),Eq (60),,' floatformat],Lminbap);
    
    fprintf(fid_log,['Lbda (dB),Eq (61),,' floatformat],Lbda(pol));
    fprintf(fid_log,['Lbam (dB),Eq (62),,' floatformat],Lbam(pol));
    fprintf(fid_log,['Lbs (dB),Eq (44),,' floatformat],Lbs);
    %     fprintf(fid_log,['Lbu (dB);Eq (63);;' floatformat],Lbu);
    %     fprintf(fid_log,['Aht (dB);Eq (64);;' floatformat],Aht);
    %     fprintf(fid_log,['Ahr (dB);Eq (64);;' floatformat],Ahr);
    fprintf(fid_log,['Lbc (dB),Eq (63),,' floatformat],Lbc);
    fprintf(fid_log,['Lb (dB),Eq (69),,' floatformat],Lb);
    fprintf(fid_log,['Ep (dBuV/m),Eq (70),,' floatformat],Ep);
    fprintf(fid_log,['Ep (dBuV/m) w.r.t. Ptx,,,' floatformat],EpPtx);
end

Ep = EpPtx;

if (debug==1)
    
    if(inside_file)
        try
            fclose(fid_log);
        end
    end
    
end

return
end