function [E,K] = energy_from_flag(all_m,empirical_topo,flag,D)
% [E,K] = energy_from_flag(all_m,empirical_topo,flag,D)
%  Compute the energy of a FLaG landscape against a population of empirical
%  connectoems
%
%   input: 
%          all_m - (n_pop,1) matrix of empirical m (number of connections)
%          empirical_topo - (4,n_pop) cell of empirical topo, output from
%               "entwork_topo" function.
%          flag - 2 dimensional matrix with (m,n_lsnets), output from
%               "generate_connections" function by specifying m as the largest
%                number in the population.
%          D - group average distance matrix
%   output: 
%           E - 2 dimensional energy matrix, (n_lsnets,n_pop)
%           K - 3 dimensional KS matrix, (n_lsnets,4,n_pop)

n = length(D);
unique_m = unique(all_m);
max_m = max(unique_m);
m_flag = size(flag,1);
if max_m>m_flag
    error('network size in flag %d is smaller than network size of empirical %d',m_flag,m_max);
end

n_lsnets = size(flag,2);
n_pop = length(all_m);

K = zeros(n_lsnets,4,n_pop);

for i = 1:length(unique_m)
    temp_m = unique_m(i);
    temp_flag = flag(1:temp_m,:);

    sbj_indx = find(all_m == temp_m);
    sbj_topo = empirical_topo(:,sbj_indx);

    temp_flag_nets = zeros(n,n,n_lsnets,"logical");
    for j = 1:n_lsnets
        temp_flag_nets_j = zeros(n,n,"logical");
        temp_flag_nets_j(temp_flag(:,j)) = 1;
        temp_flag_nets_j = temp_flag_nets_j + temp_flag_nets_j';
        temp_flag_nets(:,:,j) = temp_flag_nets_j;
    end
    ls_topo = network_topo(temp_flag_nets,D);

    for ii = 1:length(sbj_indx)
        sbj = sbj_indx(ii);
        for jj = 1:n_lsnets
            for kk = 1:4
                K(jj,kk,sbj) = fcn_ks(sbj_topo{kk,ii},ls_topo{kk,jj});
            end
        end
    end
end
E = squeeze(max(K,[],2));
E = reshape(E,n_lsnets,n_pop);


end