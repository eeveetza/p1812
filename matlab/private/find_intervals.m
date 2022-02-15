function [start,stop,pathInd] = find_intervals(series)
%find_intervals Find all intervals with consecutive trues
%     [start,stop] = find_intervals(series)
%     This function finds all true intervals, namely, the indices when the
%     intervals start and where they end
%
%     For example, for the input indices
%           0 0 1 1 1 1 0 0 0 1 1 0 0
%           0 1 1 1 1 1 1 0 1 1 1 1 0
%     this function will give back
%       start = [3;10;2;9]
%       stop = [6;11;7;12]
%       pathInd = [1;1;2;2];
%
%     Input arguments:
%     indices - Logical MxN matrix, where M is the number of paths assessed
%               in parallel and N is the number of profile points per path
%
%     Output arguments:
%     start   - vector of start-indices of the found intervals
%     stop    - vector of end-indices of the found intervals
%     pathInd - vector of path-indices of the found intervals
%
%     Example:
%     [start, stop] = find_intervals(indices)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    22JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab
%     v1    22JUN16     Roger LeClair, leclairtelecom   Modified to optimize for speed
%     v2    08MAR21     Kostas Konstantinou, Ofcom      Input 'series' can be a matrix

M = size(series,1);  % Number of paths assessed in parallel

if all(series(:)==1)
    start = ones(M,1);
    stop = size(series,2) * ones(M,1);
    pathInd = (1:M).';
else
    [start,pathInd] = find(diff([false(M,1) series],1,2).'==1);
    [stop,~] = find(diff([series false(M,1)],1,2).'==-1);
end

return
end
