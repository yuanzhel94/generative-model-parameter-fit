function [P,E,K,b] = flag_search(xrange,yrange,modeltype,modelvar,C,Aseed,D,strategy,spec,epsilon)
% function [P,E,K,flag] = flag_search(xrange,yrange,modeltype,modelvar,C,Aseed,D,strategy,spec)
%   Search parameters for generative models (spatial + topology) using FLaG
%   
%   Inputs:
%       xrange: limit of first parameter [x1,x2] (parameter eta)
%       yrange: limit of second parameter [y1,y2] (parameter gamma)
%       modeltype: generative rules as in function "generate connections"
%       modelvar: power-law or exponential-law, same as in function "generate connections"
%       C: n-by-n-by-n_pop matrix of binary unweighted empirical connectomes
%       Aseed: the seed network to initiate generative model
%       D: the group average distance matrix used for the population
%       strategy: specify which search strategy to use. Two options available. 
%           At this moment, only uniform grid search supported, use value 'grid'
%       spec: detailed specification of the search strategy as a cell array.
%           If strategy=='grid', spec = {n1,n2}, where n1 and n2 are
%               integers (how many points to sample on each parameter axis)
%           
%   Outputs:
%       P: n_points-by-2 matrix, each row is a selected point. Ordered by
%           selection sequnce if using 'vt'.
%       E: n_points-by-n_pop matrix contains the energy of each point
%           corresponds to each subject.
%       K: n_points-by-4-by-n_pop contains the four ks statistics of each
%           point corresponds to each subject.
%       b: max_m-by-n_points matrix contains the FLaG generated network with
%           maximum network density of each parameter point.
if ~exist('epsilon','var')
    epsilon = 1e-5;
end


all_m = squeeze(sum(C,[1,2]))/2;
max_m = max(all_m);
empirical_topo = network_topo(C,D);

if strcmp(strategy,'grid')
    n1 = spec{1};
    n2 = spec{2};
    eta = linspace(xrange(1),xrange(2),n1);
    gam = linspace(yrange(1),yrange(2),n2);
    eta_mat = repmat(eta,n2,1);
    gam_mat = repmat(gam',1,n1);
    P(:,1) = eta_mat(:);
    P(:,2) = gam_mat(:);

    b = generate_connections(Aseed,D,max_m,modeltype,modelvar,P,epsilon);
    [E,K] = energy_from_flag(all_m,empirical_topo,b,D);

else
    error('check "strategy" used, should be either "grid", or "vt"');
end
end


function E_ref = group_energy(Ei,rule)
if strcmp(rule,'unif')
    max_e = max(Ei,[],2);
    min_e = min(Ei,[],2);
    E_ref = min_e + (max_e - min_e).* rand(size(max_e));
elseif strcmp(rule,'mean')
    E_ref = mean(Ei,2);
elseif strcmp(rule,'median')
    E_ref = median(Ei,2);
end
end
