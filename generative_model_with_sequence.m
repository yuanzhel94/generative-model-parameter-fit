function b = generative_model_with_sequence(A,D,m,params,epsilon)
% matching index model, connections ordered by formation sequence (to enable dependent network strategy)
% warning: dependency on Brain Connectivity Toolbox
% modified based on the work Betzel et al (2016)
% "Generative models of the human connectome" in NeuroImage
% Inputs:
%         A,        binary network of seed connections
%         D,        group average Euclidean distane matrix
%         m,        largest number of network connections present in the evaluated cohort
%         params,   parameters evaluated, each row is a parameter combination
%         epsilon,  baseline probability of forming connections, to avoid undefined probability. default = 1e-5
% Output:
%         b, m x number of networks matrix of connections. 
%            m connections are order by [seed connections; generated connections by formation sequence].
    if ~exist('epsilon','var')
        epsilon = 1e-5;
    end

    n = length(D);
    nparams = size(params,1);
    b = zeros(m,nparams);

    Kseed = matching_ind(A);
    Kseed = Kseed + Kseed';
    for iparam = 1:nparams
        eta = params(iparam,1);
        gam = params(iparam,2);
        b(:,iparam) = fcn_matching(A,Kseed,D,m,eta,gam,epsilon);
    end

end

function b = fcn_matching(A,K,D,m,eta,gam,epsilon)
    b = zeros(m,1);
    K = K + epsilon;
    n = length(D);
    mseed = nnz(A)/2;
    if mseed ~= 0 %modified to record seed connections at the begining of b
        b(1:mseed) = find(triu(A,1));
    end
    Fd = D.^eta;
    Fk = K.^gam;
    Ff = Fd.*Fk.*~A;
    [u,v] = find(triu(ones(n),1));
    indx = (v - 1)*n + u;
    P = Ff(indx);
    for ii = (mseed + 1):m
        C = [0; cumsum(P)];
        r = sum(rand*C(end) >= C);
        uu = u(r);
        vv = v(r);

        A(uu,vv) = 1;
        A(vv,uu) = 1;
        b(ii) = (vv - 1)*n + uu; %modified to update b after each generative step, ensure ordered by formation sequence
        
        
        updateuu = find(A*A(:,uu));
        updateuu(updateuu == uu) = [];
        updateuu(updateuu == vv) = [];

        updatevv = find(A*A(:,vv));
        updatevv(updatevv == uu) = [];
        updatevv(updatevv == vv) = [];

        c1 = [A(:,uu)', A(uu,:)];
        for i = 1:length(updateuu)
            j = updateuu(i);
            c2 = [A(:,j)' A(j,:)];
            use = ~(~c1&~c2);
            use(uu) = 0;  use(uu+n) = 0;
            use(j) = 0;  use(j+n) = 0;
            ncon = sum(c1(use))+sum(c2(use));
            if (ncon==0)
                K(uu,j) = epsilon;
                K(j,uu) = epsilon;
            else
                K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                K(j,uu) = K(uu,j);
            end

        end

        c1 = [A(:,vv)', A(vv,:)];
        for i = 1:length(updatevv)
            j = updatevv(i);
            c2 = [A(:,j)' A(j,:)];
            use = ~(~c1&~c2);
            use(vv) = 0;  use(vv+n) = 0;
            use(j) = 0;  use(j+n) = 0;
            ncon = sum(c1(use))+sum(c2(use));
            if (ncon==0)
                K(vv,j) = epsilon;
                K(j,vv) = epsilon;
            else
                K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
                K(j,vv) = K(vv,j);
            end
        end
        Fk = K.^gam;
        Ff = Fd.*Fk.*~A;
        P = Ff(indx);
    end
end
