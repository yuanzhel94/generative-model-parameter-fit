tic;
%warning: dependency on Brain Connectivity Toolbox
% this is an toy example demonstrating the use of dependent network strategy
% generate 10x10 energy landscapes for 2 subjects using dependent network strategy
n = 84; %84x84 Desikan Killiany atlas
s = zeros(84); %initiate generative algorithm without seed connections
n_var = 10; % number of parameter sampling along each parameter dimension, n_var = 100 used for the present work
% eta = linspace(-7,2,n_var);
eta = linspace(-3.6,0.1,n_var);
% gam = linspace(0,1,n_var);
gam = linspace(0.35,0.6,n_var);
load('flag_data.mat','c'); % empirical networks for 2 subjects, number of connections m1 = 326, m2 = 355.
load('flag_data.mat','D'); % inter-regional distance matrix of Desikan Killiany atlas
m = 355; %max(m1,m2)
modeltype = "matching";
modelvar = {"powerlaw","powerlaw"};


networks = zeros(m,n_var,n_var);
for i=1:n_var
    this_eta = eta(i);
    for j=1:n_var
        this_gam = gam(j);
        params = [this_eta,this_gam];
        b = flag(s,D,m,modeltype,modelvar,params);
        networks(:,j,i) = b;
    end
end
%optionally save generated synthetic networks
% save('networks.mat','networks');

Energy = zeros(n_var,n_var,size(c,3));

for sbj = 1:size(c,3) %loop over all individuals
    this_c = squeeze(c(:,:,sbj)); %corresponding connectome
    this_m = nnz(this_c)/2; %corresponding number of connections
    this_networks = networks(1:this_m,:,:); %re-build networks with corresponding number of connections
    
    %compute topological distribution of this empirical network
    x = cell(4,1);
    x{1} = sum(this_c,2);
    x{2} = clustering_coef_bu(this_c);
    x{3} = betweenness_bin(this_c)';
    x{4} = D(triu(this_c,1) > 0);

    %compute topological distribution of each synthetic network and corresponding energy
    for i = 1:n_var
        for j=1:n_var
            b = this_networks(:,j,i);
            B = zeros(n);
            B(b) = 1;
            B = B + B';

            y = cell(4,1);
            y{1} = sum(B,2);
            y{2} = clustering_coef_bu(B);
            y{3} = betweenness_bin(B)';
            y{4} = D(triu(B,1) > 0);

            K = zeros(4,1);
            for k = 1:4
                K(k) = fcn_ks(x{k},y{k});
            end
            Energy(j,i,sbj) = max(K);
        end
    end
end

%visualize the energy landscape for selected subjects
figure;imshow(squeeze(Energy(:,:,1)),'InitialMagnification',1600);colormap(jet);colorbar; axis square;
figure;imshow(squeeze(Energy(:,:,2)),'InitialMagnification',1600);colormap(jet);colorbar; axis square;
%clearly see dependency between energy landscapes of the two individuals using FLaG.
%this is because the optimal parameters the resolution of 10x10 parameter landscape is too coarse
%but we show its effects on parameter estimation is limited when the resolution of energy
%landscapes are large (more sampling points) and using the multi-landscapes method (large K)
toc;

function kstat = fcn_ks(x1,x2) %adapted from Betzel et al (2016).
binEdges    =  [-inf ; sort([x1;x2]) ; inf];

binCounts1  =  histc (x1 , binEdges, 1);
binCounts2  =  histc (x2 , binEdges, 1);

sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);

sampleCDF1  =  sumCounts1(1:end-1);
sampleCDF2  =  sumCounts2(1:end-1);

deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
kstat = max(deltaCDF);
end