function   omega = path_fraction(d, zone, zone_r)
%path_fraction Path fraction belonging to a given zone_r
%     omega = path_fraction(d, zone, zone_r)
%     This function computes the path fraction belonging to a given zone_r
%     of the great-circle path (km) 
%
%     Input arguments:
%     d       -   matrix of distances in the path profile
%     zone    -   matrix of zones in the path profile
%     zone_r  -   reference zone for which the fraction is computed
%
%     Output arguments:
%     omega   -   path fraction belonging to the given zone_r
%
%     Example:
%     omega = path_fraction(d, zone, zone_r)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    02FEB16     Ivica Stevanovic, OFCOM         First implementation in matlab
%     v1    08MAR21     Kostas Konstantinou, Ofcom      d and zone can be MxN matrices, where M is the number of paths and N is the number of points per path

[start,stop,pathInd] = find_intervals(zone==zone_r);

if isempty(start)
    dm = zeros(size(d,1),1,class(d));
else
    delta = zeros(size(start),class(d));
    if size(d,1) > 1
        IND = d(sub2ind(size(d),pathInd,stop)) < d(pathInd,end);
    else
        IND = d(stop).' < d(pathInd,end);
    end
    if any(IND)
        ind1 = sub2ind(size(d),pathInd(IND),stop(IND)+1);
        ind2 = sub2ind(size(d),pathInd(IND),stop(IND));
        if size(d,1) > 1
            delta(IND) = delta(IND) + (d(ind1)-d(ind2))./2.0;
        else
            delta(IND) = delta(IND) + (d(ind1)-d(ind2)).'./2.0;
        end
    end
    IND = d(sub2ind(size(d),pathInd,start)) > 0;
    if any(IND)
        ind1 = sub2ind(size(d),pathInd(IND),stop(IND));
        ind2 = sub2ind(size(d),pathInd(IND),stop(IND)-1);
        if size(d,1) > 1
            delta(IND) = delta(IND) + (d(ind1)-d(ind2))./2.0;
        else
            delta(IND) = delta(IND) + (d(ind1)-d(ind2)).'./2.0;
        end
    end
    ind1 = sub2ind(size(d),pathInd,stop);
    ind2 = sub2ind(size(d),pathInd,start);
    if size(d,1) > 1
        dm = accumarray(pathInd,d(ind1)-d(ind2)+delta,[size(d,1) 1]);
    else
        dm = sum(d(ind1).'-d(ind2).'+delta);
    end
end

omega = dm ./ (d(:,end)-d(:,1));

return
end