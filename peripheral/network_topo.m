function topo = network_topo(nets,D)
% function topo = network_topo(nets,D)
%   This function compute four topological distributions of networks
%   The four topological distributions are:
%       1) degree distribution (k)
%       2) clustering coefficient distribution (c)
%       3) edge betweenness centrality distribution (e)
%       4) Euclidean space edge length distribution (d)
%
%   Inputs:
%       nets: (n,n,n_pop) binary undirected connectivity matrix, where n is
%           the number of nodes, n_pop is number of networks (subjects)
%       D: the Euclidean space distance matrix used for the population
%
%   Output:
%       topo: (n_pop,4) cell, where each row contains the four topological
%           distributions for a subject (ordered as k, c, e, d).

n_pop = size(nets,3);
topo = cell(4,n_pop);
for i = 1:n_pop
    net = squeeze(nets(:,:,i));
    topo{1,i} = sum(net,2);
    topo{2,i} = clustering_coef_bu(net);
    topo{3,i} = betweenness_bin(net)';
    topo{4,i} = D(triu(net,1) > 0);
end

end