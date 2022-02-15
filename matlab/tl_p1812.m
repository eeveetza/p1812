function [Lb, Ep] = tl_p1812(f, p, d, h, R, Ct, zone, htg, hrg, phi_t, phi_r, lam_t, lam_r, pol, varargin)
%tl_p1812 basic transmission loss according to P.1812-6
%   [Lb Ep] = tl_p1812(f, p, d, h, zone, htg, hrg, phi_t, phi_r, lam_t, lam_r, pol, varargin)
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
%     phi_t   -   Latitude of Tx station (degrees)
%     phi_r   -   Latitude of Rx station (degrees)
%     lam_t   -   Longitude of Tx station (degrees)
%     lam_r   -   Longitude of Rx station (degrees)
%     pol     -   polarization of the signal (1) horizontal, (2) vertical
%     pL      -   Required time percentage for which the calculated basic
%                 transmission loss is not exceeded (1% - 99%)
%     sigmaLoc-   location variability standard deviations computed using
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
%     Output parameters:
%     Lb     -   basic  transmission loss according to P.1812-6
%     Ep     -   the field strength relative to Ptx
%
%     Example:
%     Lb = tl_p1812(f, p, d, h, R, Ct, zone, htg, hrg, phi_t, phi_r, pol, dct, dcr, DN, N0)
%     

% How to use:
% 
% 1) by invoking only the first fourteen required input arguments:
% 
%   [Lb, Ep] = tl_p1812(f, p, d, h, R, Ct, zone, htg, hrg, phi_t, phi_r, lam_t, lam_r, pol)
%
% 2) by explicitly invoking all the input arguments:
%
%   [Lb, Ep] = tl_p1812(f, p, d, h, R, Ct, zone, htg, hrg, phi_t, phi_r, lam_t, lam_r, pol, ...
%                       pL, sigmaLoc, Ptx, DN, N0, dct, dcr, flag4, debug, fid_log);
%
% 3) or by explicitly omitting some of the optional input arguments:
% 
%   [Lb, Ep] = tl_p1812(f, p, d, h, R, Ct, zone, htg, hrg, phi_t, phi_r, lam_t, lam_r, pol, ...
%                       pL, sigmaLoc, [], DN, N0, [], [], flag4, debug, fid_log);
%
% Numbers refer to Rec. ITU-R P.1812-6

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    02JAN16     Ivica Stevanovic, OFCOM         Initial version for P.452
%     v1    07JUL16     Ivica Stevanovic, OFCOM         Modified for P.1812 
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
%     v9    11FEB22     Ivica Stevanovic, OFCOM         Aligned with P.1812-6, renamed subfolder "src" into "private" 
%                                                       which is automatically in MATLAB search path ..
%                                                       
%     

% MATLAB Version 9.2.0.556344 (R2017a) used in development of this code
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
% IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
% OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% OTHER DEALINGS IN THE SOFTWARE.
%
% THE AUTHOR(S) AND OFCOM (CH) DO NOT PROVIDE ANY SUPPORT FOR THIS SOFTWARE
%
% This function calls other functions that are placed in the ./private folder


%% Read the input arguments and check them
Nreqmax = 24;
Nreqmin = 14;

if nargin > Nreqmax    warning(['Too many input arguments; The function requires at most ' num2str(Nreqmax) ,...
        'input arguments. Additional values ignored. Input values may be wrongly assigned.']);
end

if nargin < Nreqmin
    error(['Function requires at least ' num2str(Nreqmin) ' input arguments.']);
end
%

% TODO: make sure the necessary subroutines are in the Octave path
%     s = pwd;
%     if ~exist('dl_bull.m','file')
%         addpath([s '/src/'])
%     end
 

    
    % verify input argument values and limits
    check_limit(f, 0.03, 6.0, 'f [GHz]');
    check_limit(p, 1, 50, 'p [%]');
    check_limit(phi_t, -80, 80, 'phi_t [deg]');
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

% Optional arguments (default values)
pL = 50;
sigma = 0;
Ptx = 1;
DN = 45; 
N0 = 325;
dct = 500;
dcr = 500;
flag4 = 0;

if zone(1) == 1 % Tx at sea
 dct = 0;
end
if zone(end) == 1 %Rx at sea
 dcr = 0;
end

% ws = 27; 
debug = 0;
fid_log = [];

icount = Nreqmin+ 1;

if nargin >=icount
    if (~isempty(varargin{1}))
        pL=varargin{1};
    end
    if nargin >=icount+1
        if (~isempty(varargin{2}))
            sigmaLoc=varargin{2};
        end
        if nargin >=icount+2
            if (~isempty(varargin{3}))
                Ptx=varargin{3};
            end
            if nargin >=icount+3
                if(~isempty(varargin{4}))
                    DN=varargin{4};
                end
                if nargin >=icount + 4
                    if (~isempty(varargin{5}))
                        N0=varargin{5};
                    end
                    if nargin >=icount + 5
                        if (~isempty(varargin{6}))
                            dct=varargin{6};
                        end
                        if nargin >=icount + 6
                            if(~isempty(varargin{7}))
                                dcr=varargin{7};
                            end
                            if nargin >=icount + 7
                                if(~isempty(varargin{8}))
                                    flag4=varargin{8};
                                end
                                if nargin >=icount + 8
                                    if(~isempty(varargin{9}))
                                        debug=varargin{9};
                                    end
                                    if nargin >=icount + 9
                                        if(~isempty(varargin{10}))
                                            fid_log=varargin{10};
                                        end
                                        
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
        end
    end
end
check_limit(pL, 1, 99, 'pL [%]');

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

floatformat= '%.10g;\n';

if (debug)
    
    fprintf(fid_log,'# Parameter;Ref;;Value;\n');
    fprintf(fid_log,['Ptx (kW);;;' floatformat],Ptx);
    fprintf(fid_log,['f (GHz);;;' floatformat],f);
    fprintf(fid_log,['p (%%);;;' floatformat],p);
    fprintf(fid_log,['pL (%%);;;' floatformat],pL);
    fprintf(fid_log,['sigma (dB);;;' floatformat],sigma);
    fprintf(fid_log,['phi_t (deg);;;' floatformat],phi_t);
    fprintf(fid_log,['phi_r (deg);;;' floatformat],phi_r);
    fprintf(fid_log,['lam_t (deg);;;' floatformat],lam_t);
    fprintf(fid_log,['lam_r (deg);;;' floatformat],lam_r);
    fprintf(fid_log,['htg (m);;;' floatformat],htg);
    fprintf(fid_log,['hrg (m);;;' floatformat],hrg);
    fprintf(fid_log,['pol;;;' '%d\n'],pol);
%    fprintf(fid_log,['ws (m);;;' floatformat],ws);
    fprintf(fid_log,['DN ;;;' floatformat],DN);
    fprintf(fid_log,['N0 ;;;' floatformat],N0);
    fprintf(fid_log,['dct (km) ;;;' floatformat],dct);
    fprintf(fid_log,['dcr (km) ;;;' floatformat],dcr);
    fprintf(fid_log,['R2 (m) ;;;' floatformat],R(2));
    fprintf(fid_log,['Rn-1 (m) ;;;' floatformat],R(end-1));
    fprintf(fid_log,['Ct Tx  ;Table 2;;' floatformat],Ct(2));
    fprintf(fid_log,['Ct Rx ;Table 2;;' floatformat],Ct(end-1));
    
end


% Compute the path profile parameters
% Path center latitude
% phi_path = (phi_t + phi_r)/2;
%[phi_path, ~] = gcintermediate(phi_t, lam_t, phi_r, lam_r, 0.5);
Re = 6371;
dpnt = 0.5*(d(end)-d(1));
[~, phi_path,~] = great_circle_path(lam_r, lam_t, phi_r, phi_t, Re, dpnt);

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
    fprintf(fid_log,';;;;\n');
    %fprintf(fid_log,'# Path;Ref;;Value;\n');
    fprintf(fid_log,['d (km);;;' floatformat],dtot);
    fprintf(fid_log,['dlt (km);Eq (80);;' floatformat],dlt);
    fprintf(fid_log,['dlr (km);Eq (83);;' floatformat],dlr);
    fprintf(fid_log,['th_t (mrad);Eqs (78-79);;' floatformat],theta_t);
    fprintf(fid_log,['th_r (mrad);Eqs (81-82);;' floatformat],theta_r);
    fprintf(fid_log,['th (mrad);Eq (84);;' floatformat],theta);
    fprintf(fid_log,['hts (m);;;' floatformat],hts);
    fprintf(fid_log,['hrs (m);;;' floatformat],hrs);
    fprintf(fid_log,['htc (m);Table 5;;' floatformat],htc);
    fprintf(fid_log,['hrc (m);Table 5;;' floatformat],hrc);
    fprintf(fid_log,['w;Table 5;;' floatformat],omega);
    fprintf(fid_log,['dtm (km);Sec 3.6;;' floatformat],dtm);
    fprintf(fid_log,['dlm (km);Sec 3.6;;' floatformat],dlm);
    fprintf(fid_log,['phi (deg);Eq (4);;' floatformat],phi_path);
    fprintf(fid_log,['b0 (%%);Eq (5);;' floatformat],b0);
    fprintf(fid_log,['ae (km);Eq (7a);;' floatformat],ae);
    fprintf(fid_log,['hst (m);Eq (87);;' floatformat],hst_n);
    fprintf(fid_log,['hsr (m);Eq (88);;' floatformat],hsr_n);
    fprintf(fid_log,['hst (m);Eq (92a);;' floatformat],hst);
    fprintf(fid_log,['hsr (m);Eq (92b);;' floatformat],hsr);
    fprintf(fid_log,['hstd (m);Eq (91);;' floatformat],hstd);
    fprintf(fid_log,['hsrd (m);Eq (91);;' floatformat],hsrd);
    fprintf(fid_log,['htc'' (m);Eq (37a);;' floatformat],htc-hstd);
    fprintf(fid_log,['hrc'' (m);Eq (37b);;' floatformat],hrc-hsrd);
    fprintf(fid_log,['hte (m);Eq (94);;' floatformat],hte);
    fprintf(fid_log,['hre (m);Eq (94);;' floatformat],hre);
    fprintf(fid_log,['hm (m);Eq (95);;' floatformat],hm);
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

if Lminbap <= Lbd
   Lbda = Lminbap + (Lbd-Lminbap)*Fk; % eq (61)
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


Lloc = 0; % outdoors only (67a)

% Location variability of losses (Section 4.8)
if zone(end) ~= 1 %Rx at sea
    Lloc = -inv_cum_norm( pL/100 ) * sigmaLoc;
end

% Basic transmission loss not exceeded for p% time and pL% locations
% (Sections 4.8 and 4.9) not implemented

Lb = max(Lb0p, Lbc + Lloc);  %  eq (69)

% The field strength exceeded for p% time and pL% locations

Ep = 199.36+ 20*log10(f) - Lb; % eq (70)

% % Scale to the transmitter power

EpPtx = Ep + 10*log10(Ptx);

 if (debug)  
    fprintf(fid_log,['Fi;Eq (40a);;' floatformat],Fi);     
    fprintf(fid_log,['Fj;Eq (57);;' floatformat],Fj);
    fprintf(fid_log,['Fk;Eq (58);;' floatformat],Fk);
    fprintf(fid_log,['Lbfs;Eq (8);;' floatformat],Lbfs);
    fprintf(fid_log,['Lb0p;Eq (10);;' floatformat],Lb0p);
    fprintf(fid_log,['Lb0b;Eq (11);;' floatformat],Lb0b);
    fprintf(fid_log,['Lbulla (dB);Eq (21);;' floatformat],Lbulla50);
    fprintf(fid_log,['Lbulls (dB);Eq (21);;' floatformat],Lbulls50);
    fprintf(fid_log,['Ldsph (dB);Eq (27);;' floatformat],Ldsph50(pol));
    fprintf(fid_log,['Ld50 (dB);Eq (39);;' floatformat],Ld50(pol));
    fprintf(fid_log,['Ldb (dB);Eq (39);;' floatformat],Ldb(pol));    
    fprintf(fid_log,['Ldp (dB);Eq (41);;' floatformat],Ldp(pol));
    fprintf(fid_log,['Lbd50 (dB);Eq (42);;' floatformat],Lbd50(pol));
    fprintf(fid_log,['Lbd (dB);Eq (43);;' floatformat],Lbd(pol));   
    
    fprintf(fid_log,['Lminb0p (dB);Eq (59);;' floatformat],Lminb0p(pol));
    
    fprintf(fid_log,['Lba (dB);Eq (46);;' floatformat],Lba);
    fprintf(fid_log,['Lminbap (dB);Eq (60);;' floatformat],Lminbap);
    
    fprintf(fid_log,['Lbda (dB);Eq (61);;' floatformat],Lbda(pol));    
    fprintf(fid_log,['Lbam (dB);Eq (62);;' floatformat],Lbam(pol));    
    fprintf(fid_log,['Lbs (dB);Eq (44);;' floatformat],Lbs);
%     fprintf(fid_log,['Lbu (dB);Eq (63);;' floatformat],Lbu);
%     fprintf(fid_log,['Aht (dB);Eq (64);;' floatformat],Aht);
%     fprintf(fid_log,['Ahr (dB);Eq (64);;' floatformat],Ahr);
    fprintf(fid_log,['Lbc (dB);Eq (65);;' floatformat],Lbc);
    fprintf(fid_log,['Lb (dB);Eq (71);;' floatformat],Lb);
    fprintf(fid_log,['Ep (dBuV/m);Eq (72);;' floatformat],Ep);
    fprintf(fid_log,['Ep (dBuV/m) w.r.t. Ptx;;;' floatformat],EpPtx);
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