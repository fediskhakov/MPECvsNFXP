function [logl, ev, score]=ll(data, X, P, mp, par, method)
% * llu: Partial likelihood function, discrete choice
% Syntax: [ll, ev, score]=llu(data, nod, P, mp, par)
mp.p=abs(mp.p); %helps BHHH which is run as unconstrained optimization
n_c=numel(mp.c);
N=size(data.x,1);

% Update u, du and P evaluated in grid points
dc=0.001*X.grid;
cost=mp.c*0.001*X.grid;
if strcmp(method,'fullmle')
    P = nfxp.statetransition(mp.p, X.n);
end

% Solve model
[ev, pk, F]=nfxp.solve(mp.ev, P, cost, mp, par);
% Evaluate likelihood function
% log likelihood regarding replacement choice
lp=pk(data.x);
logl=log(lp+(1-2*lp).*(data.d==2));

% fprintf('ll params [%1.5f %1.5f %1.5f %1.5f %1.5f %1.5f]\n',mp.c,mp.p,mp.RC);
% fprintf('ll min pk=%1.20f max pk=%1.20f\n',min(pk),max(pk))
% plot (pk)
% hold on
% pause

% add on log like for mileage process
if strcmp(method,'fullmle')
    p=[mp.p; 1-sum(mp.p)];
    n_p=numel(p)-1;
    logl=logl + log(p(1+ data.dx1));
else
    n_p=0;
end


if nargout ==3;  %% compute scores
    
    % step 1: compute derivative of contraction operator wrt. parameters
    dtdmp=zeros(X.n,1+ n_c + n_p);
    dtdmp(:, 1)=-P*(1-pk);
    dtdmp(:, 2:1+n_c)=-(P*pk).*dc; %
    if strcmp(method,'fullmle')
        dtp=log(exp(-cost+mp.beta*mp.ev)+exp(-mp.RC-cost(1)+mp.beta*mp.ev(1)));
        for iP=1:n_p;
            dtdmp(1:X.n-iP, 1+n_c+iP)= [dtp(iP:end-1)] - [dtp(n_p+1:X.n); repmat(dtp(end), n_p-iP, 1)];
        end
        invp=exp(-log(p));
        invp=[sparse(1:n_p,1:n_p,invp(1:n_p),n_p,n_p); -ones(1,n_p)*invp(n_p+1)];
        N=size(data.x,1);
    end
    
    % step 2: compute derivative of ev wrt. parameters
    devdmp=F\dtdmp;  
    
    % step 3: compute derivative of log-likelihood wrt. parameters
    score=bsxfun(@times, (lp- 1 + (data.d==2)),[-ones(N,1) dc(data.x,:) zeros(N,n_p)] + (devdmp(ones(N,1),:)-devdmp(data.x,:)));    
    if strcmp(method,'fullmle')
        for iP=1:n_p;
            score(:,1+n_c+iP)= score(:,1+n_c+iP) + invp(1+data.dx1,iP);
        end
    end
end


end % end of ll