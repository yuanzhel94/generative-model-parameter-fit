function [est,avgE] = multilandscape_estimate(P,E)
%function est = multilandscape_estimate(P,E)
%   This function compute the multilandscape estimate for each subject from
%   the sampled points from each landscape, and their energies.
%       Inputs:
%           P: n_points-by-2-by-K_ls matrix, where K_ls is the number of
    %           landscapes will be evaluated. P(i,:,j) contains the sampled
    %           parameter points for the i-th point sampled in landscape j.
%           E: n_points-by-n_pop-by-K_ls matrix. E(i,j,k) contains the
%               energy of the i-th point in landscape k against empirical
%               network j.

%       Output:
%           est: n_pop-by-2 matrix, where each row contains the estimated
%               model parameters of a subject.
%           avgE: n_pop-by-1 vector, each element is the average energy for
%               all points selected for the estimation.

K_ls = size(P,3);
n_pop = size(E,2);
temp_est = zeros(n_pop,2,K_ls);
Emin = zeros(n_pop,K_ls);
for k = 1:K_ls
    Ek = squeeze(E(:,:,k));
    for sbj = 1:n_pop
        Emin(sbj,k) = min(Ek(:,sbj));
        indx = find(Ek(:,sbj)==Emin(sbj,k));
        temp_est(sbj,:,k) = mean(P(indx,:,k),1);
    end
end
est = mean(temp_est,3);
avgE = mean(Emin,2);
end