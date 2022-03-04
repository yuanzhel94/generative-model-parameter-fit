% This script estimates the sample sizes required for an expected
% between-group differences in eta or gamma. Or do the reverse estimation.
% This is achieved using the individual variation in estiamtes for the HCP-YA cohort (std)
% Use limited to connectome mapped to the Desikan Killiany atals, 10% network density on average, and generative algorithms initiated without seed connections.
load('sample_size_estimate_data.mat');
Power = 0.8;
T=norminv(1-0.025);
k_indx = 3; %select K, 1: K = 1; 2: K = 10; 3: K = 50; 4: K = 100.
eta_sigma1 = eta_std(k_indx);eta_sigma = sqrt(2*eta_sigma1^2);
gam_sigma1 = gam_std(k_indx);gam_sigma = sqrt(2*gam_sigma1^2);
%Either run part 1 or part 2 for different purpose


%PART 1
%estimate minimum delta resolvable with a given sample size
samp = 200; %given sample size

%minimum delta in eta
eta_delta = (T-norminv(1-Power)) * eta_sigma /sqrt(samp);
%minimum delta in gamma
gam_delta = (T-norminv(1-Power)) * gam_sigma /sqrt(samp);


%PART2
%estimate minimum sample size to resolve a given delta in eta or gamma
eta_delta = 0.028;
samp1 = ((T-norminv(1-Power))*eta_sigma/eta_delta)^2;
gam_delta = 0.0037;
samp2 = ((T-norminv(1-Power))*gam_sigma/gam_delta)^2;
