function   dm = longest_cont_dist(d, zone, zone_r)
%longest_cont_dist Longest continuous path belonging to the zone_r
%     dm = longest_cont_dist(d, zone, zone_r)
%     This function computes the longest continuous section of the
%     great-circle path (km) for a given zone_r
%
%     Input arguments:
%     d       -   MxN matrix of distances in the path profile
%     zone    -   MxN matrix of zones in the path profile
%                 4 - Inland, 3 - Coastal, 1 - See
%     zone_r  -   reference zone for which the longest continuous section
%                 is computed, zone_r = 34 for combined inland-coastal land
%     ,where M is the number of paths assessed in parallel and N is the
%     number of profile points per path.
%     
%     Output arguments:
%     dm      -   the longest continuous section of the great-circle path (km) for a given zone_r
%
%     Example:
%     dm = longest_cont_dist(d, zone, zone_r)

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    22JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab
%     v1    12FEB16     Ivica Stevanovic, OFCOM         included zone_r==12
%     v2    08JUL16     Ivica Stevanovic, OFCOM         modified mapping to  GlobCover data format
%                                                       before: 2 - Inland, 1 - Coastal land, 3 - Sea
%                                                       now:    4 - Inland, 3 - Coastal land, 1 - Sea
%     v3    08MAR21     Kostas Konstantinou, Ofcom      d and zone can be MxN matrices, where M is the number of paths and N is the number of points per path

if zone_r  == 34 % inland + coastal land
    [start,stop,pathInd] = find_intervals(or(zone==3,zone==4));
else
    [start,stop,pathInd] = find_intervals(zone==zone_r);
end

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
        dm = accumarray(pathInd,d(ind1)-d(ind2)+delta,[size(d,1) 1],@max);
    else
        dm = max(d(ind1).'-d(ind2).'+delta);
    end
end

return
end