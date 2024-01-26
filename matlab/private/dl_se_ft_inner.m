function Ldft = dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f)
%%dl_se_ft_inner The inner routine of the first-term spherical diffraction loss
%   This function computes the first-term part of Spherical-Earth diffraction
%   loss exceeded for p% time for antenna heights
%   as defined in Sec. 4.3.3 of the ITU-R P.1812-6, equations (29-36)
%
%     Input parameters:
%     epsr    -   Relative permittivity
%     sigma   -   Conductivity (S/m)
%     d       -   Great-circle path distance (km)
%     hte     -   Effective height of interfering antenna (m)
%     hre     -   Effective height of interfered-with antenna (m)
%     adft    -   effective Earth radius (km)
%     f       -   frequency (GHz)
%
%     Output parameters:
%     Ldft   -   The first-term spherical-Earth diffraction loss not exceeded for p% time
%                implementing equations (29-36), Ldft(1) is for horizontal
%                and Ldft(2) for the vertical polarization
%
%     Example:
%     Ldft = dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f)

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation (P.452)
%     v1    06JUL16     Ivica Stevanovic, OFCOM         Modifications for P.1812


%% 

% Normalized factor for surface admittance for horizontal (1) and vertical
% (2) polarizations

K(1) = 0.036* (adft*f).^(-1/3) * ( (epsr-1).^2 + (18*sigma/f).^2 ).^(-1/4);   % Eq (29a)

K(2) = K(1) * (epsr.^2 + (18*sigma/f).^2).^(1/2);       % Eq (29b)

% Earth ground/polarization parameter

beta_dft = (1 + 1.6*K.^2 + 0.67* K.^4)./( 1 + 4.5* K.^2 + 1.53* K.^4);  % Eq (30)

% Normalized distance

X = 21.88* beta_dft .* (f./ adft.^2).^(1/3) * d;          % Eq (31)

% Normalized transmitter and receiver heights

Yt = 0.9575* beta_dft * (f.^2 / adft).^(1/3) * hte;       % Eq (32a)

Yr = 0.9575* beta_dft * (f.^2 / adft).^(1/3) * hre;       % Eq (32b)


Fx = [0,0];
GYt= [0,0];
GYr= [0,0];
% Calculate the distance term given by:

for ii = 1:2
    if X(ii) >= 1.6
        Fx(ii) = 11 + 10*log10(X(ii)) - 17.6*X(ii);
    else
        Fx(ii) = -20*log10(X(ii)) - 5.6488* (X(ii)).^1.425;     % Eq (33)
    end
end

Bt = beta_dft  .* Yt;             % Eq (35)

Br = beta_dft .* Yr;              % Eq (35)

for ii = 1:2
    if Bt(ii)>2
        GYt(ii) = 17.6*(Bt(ii) - 1.1).^0.5 - 5*log10(Bt(ii) -1.1)-8;
    else
        GYt(ii) = 20*log10(Bt(ii) + 0.1* Bt(ii).^3);
    end
    
    if Br(ii)>2
        GYr(ii) = 17.6*(Br(ii) - 1.1).^0.5 - 5*log10(Br(ii) -1.1)-8;
    else
        GYr(ii) = 20*log10(Br(ii) + 0.1* Br(ii).^3);
    end
    
    if GYr(ii) < 2 + 20*log10(K(ii))
        GYr(ii) = 2 + 20*log10(K(ii));
    end
    
    if GYt(ii) < 2 + 20*log10(K(ii))
        GYt(ii) = 2 + 20*log10(K(ii));
    end
    
end

Ldft = -Fx - GYt - GYr;         % Eq (36)


% floatformat= '%.10g;\n';
% fid = fopen('Ldsph.csv', 'a');
% % fprintf(fid,['adft ;;;' floatformat],adft);
% % fprintf(fid,['epsr ;;;' floatformat],epsr);
% % fprintf(fid,['sigma ;;;' floatformat],sigma);
% % fprintf(fid,['f ;;;' floatformat],f);
% fprintf(fid,['Kh ;Eq (29a);;' floatformat],K(1));
% fprintf(fid,['Kv ;Eq (29b);;' floatformat],K(2));
% fprintf(fid,['beta_dft hor ;Eq (30);;' floatformat],beta_dft(1));
% fprintf(fid,['beta_dft ver ;Eq (30);;' floatformat],beta_dft(2));
% fprintf(fid,['X hor ;Eq (31);;' floatformat],X(1));
% fprintf(fid,['X ver ;Eq (31);;' floatformat],X(2));
% fprintf(fid,['Yt hor ;Eq (32a);;' floatformat],Yt(1));
% fprintf(fid,['Yt ver ;Eq (32a);;' floatformat],Yt(2));
% fprintf(fid,['Yr hor ;Eq (32b);;' floatformat],Yr(1));
% fprintf(fid,['Yr ver ;Eq (32b);;' floatformat],Yr(2));
% fprintf(fid,['Fx hor ;Eq (33);;' floatformat],Fx(1));
% fprintf(fid,['Fx ver ;Eq (33);;' floatformat],Fx(2));
% fprintf(fid,['GYt hor ;Eq (34);;' floatformat],GYt(1));
% fprintf(fid,['GYt ver ;Eq (34);;' floatformat],GYt(2));
% fprintf(fid,['GYr hor ;Eq (34);;' floatformat],GYr(1));
% fprintf(fid,['GYr ver ;Eq (34);;' floatformat],GYr(2));
% fprintf(fid,['Bt hor ;Eq (35);;' floatformat],Bt(1));
% fprintf(fid,['Bt ver ;Eq (35);;' floatformat],Bt(2));
% fprintf(fid,['Br hor ;Eq (35);;' floatformat],Br(1));
% fprintf(fid,['Br ver ;Eq (35);;' floatformat],Br(2));
% fprintf(fid,['Bt hor ;Eq (35);;' floatformat],Bt(1));
% fprintf(fid,['Bt ver ;Eq (35);;' floatformat],Bt(2));
% fprintf(fid,['Ldft hor ;Eq (36);;' floatformat],Ldft(1));
% fprintf(fid,['Ldft ver ;Eq (36);;' floatformat],Ldft(2));
% fprintf(fid,[';;;\n']);
% fclose(fid);

return
end