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
%           If use grid search, specify 'grid'. 
%           If use Voronoi Tessellation, specify 'vt'.
%       spec: detailed specification of the search strategy as a cell array.
%           If strategy=='grid', spec = {n1,n2}, where n1 and n2 are
%               integers (how many points to sample on each parameter axis)
%           If strategy=='vt', spec={n_draw,pow,rule}
%                   n_draw (size 1-by-n_iter) and pow (size 1-by-(niter-1))are vectors,
    %                   where n_iter is the number of iterations for vt selection. n_draw contains the number of
    %                   points selected at each iteration, and pow controls the
    %                   scale factor for selection probability (e.g., 5 iterations,
    %                   5 values of n_draw, 4 values of pow).
%                   rule (either 'unif',or 'mean') controls how to compute the group-representative
%                       energy (one generated network has different energeis to different empirical networks).
%                       The group-representative energy will be used to update the sampling probability in the next iteration.
%                       If 'unif', uniformly sample a energy between the
%                       maximum and minimum energy.
%                       If 'mean', compute the mean energy for the
%                       population.
%                       If 'median', compute the median energy for the
%                       population.
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
elseif strcmp(strategy,'vt')
    n_draw = spec{1};
    pow = spec{2};
    rule = spec{3};
    n_sum = sum(n_draw);
    n_pop = size(C,3);
        
    P = zeros(n_sum,2);
    b = zeros(max_m,n_sum);
    E = zeros(n_sum,n_pop);
    K = zeros(n_sum,4,n_pop);

    etai = unifrnd(xrange(1),xrange(2),n_draw(1),1);
    gami = unifrnd(yrange(1),yrange(2),n_draw(1),1);
    Pi = [etai,gami];
    bi = generate_connections(Aseed,D,max_m,modeltype,modelvar,Pi,epsilon);
    [Ei,Ki] = energy_from_flag(all_m,empirical_topo,bi,D);
    P(1:n_draw(1),:) = Pi;
    b(:,1:n_draw(1)) = bi;
    E(1:n_draw(1),:) = Ei;
    K(1:n_draw(1),:,:) = Ki;
    E_ref = group_energy(Ei,rule);
    count = n_draw(1);
    for i = 2:length(n_draw)
        Pi = fcn_voronoi_select(Pi,E_ref,n_draw(i),xrange,yrange,pow(i-1));
        bi = generate_connections(Aseed,D,max_m,modeltype,modelvar,Pi,epsilon);
        [Ei,Ki] = energy_from_flag(all_m,empirical_topo,bi,D);
        P(count+1:count+n_draw(i),:) = Pi;
        b(:,count+1:count+n_draw(i)) = bi;
        E(count+1:count+n_draw(i),:) = Ei;
        K(count+1:count+n_draw(i),:,:) = Ki;
        E_ref = group_energy(Ei,rule);
        count = count+n_draw(i);
    end
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
