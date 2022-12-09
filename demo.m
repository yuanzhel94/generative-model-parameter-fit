%% 
addpath(genpath('~/BCT/2019_03_03_BCT/')); %change to your own Brain Connectivity Toolbox (BCT) path
addpath(genpath('./peripheral/'));
% This is a demo to introduce the use of Fast Landscape Generation (FLaG)
% in "Parameter estimation for connectome generative models: Accuracy,
% reliability, and a fast parameter fitting method".
% 
% In this demo, the generative model parameters of a group of 100 networks
% will be estimated using FLaG and Multi-point method (K = 50). It may take
% up to 3 days to run the code on a single core.
% To reduce the time of testing the code, reduce K and spec.
% Alternatively, implement your own computation in parallel.
% The generative model used in this demo is the "matching index" model
% developed in Betzel et al (2016).

% load data of 100 networks to be evaluated.
% to avoid breaching the rule for data sharing, here we use 100 networks
% generated from the matching index rule (ground truth parameters are
% stored in the variable gt_params, network sizes are sampled from the
% empirical network sizes to achieve a similar mean/std to empirical data).
load('demo_data.mat');
% number of landscapes (parameter for multi-landscape method for parameter
% estimation).
K_ls = 50;
% parameter range of interest (search in this range)
eta_range = [-7,1];
gam_range = [-0.3,1];

%specify strategy and spec
%100*100 grid

strategy = 'grid';
spec = {100,100};
n_points = spec{1}*spec{2};


%initialize for storage
n_pop = size(C,3);
max_m = max(samp_m);
P = zeros(n_points,2,K_ls);
E = zeros(n_points,n_pop,K_ls);
K = zeros(n_points,4,n_pop,K_ls);
b = zeros(max_m,n_points,K_ls);


tic;
%build the landscapes
for i = 1:K_ls
    [Pi,Ei,Ki,bi] = flag_search(eta_range,gam_range,modeltype,modelvar,C,Aseed,D,strategy,spec);
    P(:,:,i) = Pi;
    E(:,:,i) = Ei;
    K(:,:,:,i) = Ki;
    b(:,:,i) = bi;
end
[est,avgE] = multilandscape_estimate(P,E);


toc;

save demo_result.mat est avgE gt_params 

% figure;
% plot(gt_params(:,1),est(:,1));