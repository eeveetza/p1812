function I = inv_cum_norm( x )
%inv_cum_norm approximation to the inverse cummulative normal distribution
%   I = inv_cum_norm( x )
%   This function implements an approximation to the inverse cummulative
%   normal distribution function for 0<x<1 as defined in Attachment 2 to
%   Annex 1 of the ITU-R P.1812-4

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v1    09DEC16     Ivica Stevanovic, OFCOM         Version P.1812-4
%     v0    02JAN16     Ivica Stevanovic, OFCOM         Initial version P.452-16

if x < 0.000001
    x = 0.000001;
end

if x > 0.999999
    x = 0.000001;
end

if x > 0.5
    x = 1 - x;           % eq (96a)
end

tx = sqrt(-2*log(x));    % eq (97a)

C0 = 2.515516698;        % eq (97c)
C1 = 0.802853;           % eq (97d)
C2 = 0.010328;           % eq (97e)
D1 = 1.432788;           % eq (97f)
D2 = 0.189269;           % eq (97g)
D3 = 0.001308;           % eq (97h)

ksi = ( (C2*tx+C1)*tx + C0 )/ ( ((D3*tx + D2)*tx + D1)*tx + 1 );  % eq (97b)

I = ksi - tx;            % eq (96a/b)

return 
end

