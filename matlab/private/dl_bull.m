function [Lbull,maxI] = dl_bull(d,g,hts,hrs,ap,f)
%dl_bull Bullington part of the diffraction loss according to P.1812-4
%   This function computes the Bullington part of the diffraction loss
%   as defined in ITU-R P.1812-4 in 4.3.1
%
%     Input parameters:
%     d       -   vector of distances di of the i-th profile point (km)
%     g       -   vector of heights hi of the i-th profile point (meters
%                 above mean sea level) + representative clutter height 
%                 Both vectors d and g contain n+1 profile points
%     hts     -   transmitter antenna height in meters above sea level (i=0)
%     hrs     -   receiver antenna height in meters above sea level (i=n)
%     ap      -   the effective earth radius in kilometers
%     f       -   frequency expressed in GHz
%
%     Output parameters:
%     Lbull   -   Bullington diffraction loss for a given path
%     maxI    -   Path index with highest diffraction parameter
%
%     Example:
%     Lbull = dl_bull(d, g, hts, hrs, ap, f)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation for P.452-16
%     v1    06JUL16     Ivica Stevanovic, OFCOM         First implementation for P.1812-4
%     v2    09MAR21     Kostas Konstantinou, Ofcom      Input d,g can be matrices, and hts,hrs,ap can be vectors. Added maxI as output.


%% Body of function

% Effective Earth curvature Ce (km^-1)
Ce = 1 ./ ap;

% Wavelength in meters
lambda = 0.3/f;

% Complete path length
dtot = d(:,end) - d(:,1);

% Find the intermediate profile point with the highest slope of the line
% from the transmitter to the point
di = d(:,2:end-1);
gi = g(:,2:end-1);
Stim = max((gi+500.*Ce.*di.*(dtot-di)-hts)./di,[],2);  % Eq (13)

% Calculate the slope of the line from transmitter to receiver assuming a
% LoS path
Str = (hrs-hts) ./ dtot;  % Eq (14)

Luc = zeros(size(Stim),class(Stim));
IND = find(Stim<Str);
maxI = zeros(size(Stim),'uint16');
if ~isempty(IND)
    % Case 1, Path is LoS
    % Find the intermediate profile point with the highest diffraction
    % parameter nu:
    di2 = di(IND,:);
    dtot2 = dtot(IND);
    [numax,maxI(IND)] = max((gi(IND,:)+500.*Ce(IND).*di2.*(dtot2-di2) - ...
        (hts(IND).*(dtot2-di2)+hrs(IND).*di2)./dtot2) .* ...
        sqrt(0.002.*dtot2./(lambda.*di2.*(dtot2-di2))),[],2);  % Eq (15)
    
    IND2 = numax > -0.78;
    if any(IND2)
        Luc(IND(IND2)) = 6.9 + ...
            20.*log10(sqrt((numax(IND2)-0.1).^2+1)+numax(IND2)-0.1);  % Eq (12), (16)
    end
end
IND = find(~(Stim<Str));
if ~isempty(IND)
    hts2 = hts(IND);
    hrs2 = hrs(IND);
    Stim2 = Stim(IND);
    
    % Path is transhorizon
    % Find the intermediate profile point with the highest slope of the
    % line from the receiver to the point
    di2 = di(IND,:);
    dtot2 = dtot(IND);
    [Srim,maxI(IND)] = max((gi(IND,:)+500.*Ce(IND).*di2.*(dtot2-di2)-hrs2)./(dtot2-di2),[],2);  % Eq (17)
    
    % Calculate the distance of the Bullington point from the transmitter:
    dbp = (hrs2-hts2+Srim.*dtot2) ./ (Stim2+Srim);  % Eq (18)
    
    % Calculate the diffraction parameter, nub, for the Bullington point
    nub = (hts2+Stim2.*dbp-(hts2.*(dtot2-dbp)+hrs2.*dbp)./dtot2) .* ...
        sqrt(0.002.*dtot2./(lambda.*dbp.*(dtot2-dbp)));  % Eq (20)
    
    % The knife-edge loss for the Bullington point is given by
    IND2 = nub > -0.78;
    if any(IND2)
        Luc(IND(IND2)) = 6.9 + ...
            20.*log10(sqrt((nub(IND2)-0.1).^2+1)+nub(IND2)-0.1);  % Eq (12), (20)
    end
end

% For Luc calculated using either (16) or (20), Bullington diffraction loss
% for the path is given by
Lbull = Luc + (1-exp(-Luc./6.0)).*(10+0.02.*dtot);  % Eq (21)

if nargout == 2
    maxI = maxI + 1;
end
return
end
