function I = inv_cum_norm( x )
%inv_cum_norm approximation to the inverse cummulative normal distribution
%   I = inv_cum_norm( x )
%   This function implements an approximation to the inverse cummulative
%   normal distribution function for 0<x<1 as defined in Attachment 2 to
%   Annex 1 of the ITU-R P.1812-6

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v2    19MAR20     Ivica Stevanovic, OFCOM         Corrected bugs having occured when x > 0.5 
%     v1    09DEC16     Ivica Stevanovic, OFCOM         Version P.1812-4
%     v0    02JAN16     Ivica Stevanovic, OFCOM         Initial version P.452-16

if x < 0.000001
    x = 0.000001;
end

if x > 0.999999
    x = 0.999999;
end

if x<= .5
    I = T(x)-C(x);  %(96a)
else
    I = -(T(1-x)-C(1-x)); %(96b)
end
   
return

    function outT = T(y)
        outT = sqrt(-2*log(y));     %(97a)
    return

    function outC = C(z)   % (97)
        C0 = 2.515516698;
        C1 = 0.802853;
        C2 = 0.010328;
        D1 = 1.432788;
        D2 = 0.189269;
        D3 = 0.001308;
        outC = (((C2*T(z)+C1)*T(z))+C0)/(((D3*T(z)+D2)*T(z)+D1)*T(z)+1);   %(97b)
    return
