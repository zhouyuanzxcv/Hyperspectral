function [A,J2,A2] = calc_adjacency_graph_cached(Y,k)
%CALC_ADJACENCY_GRAPH_CACHED Summary of this function goes here
%   Detailed explanation goes here
[N,B] = size(Y);

cached_params = [];
cached_params.Y = Y;
cached_params.k = k;
cached_file = fullfile(tempdir, 'cached_adjacency_graph.mat');
ccr = CachedCalculationResults(cached_file);
result = ccr.GetMatchedResult(cached_params);

if isempty(result)
    [A,J2,A2] = create_adjacency_graph(Y,'nn',k,1);
    result = struct('A',A,'J2',J2,'A2',A2);
    if N*B > 1e6
        ccr.AddResult(cached_params,result);
    end
else
    A = result.A;
    J2 = result.J2;
    A2 = result.A2;
end

