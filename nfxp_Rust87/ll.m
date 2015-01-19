function [logl, ev, score, H]=ll(data, X, P, mp, par, method)
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

% add on log like for mileage process
if strcmp(method,'fullmle')
    p=[mp.p; 1-sum(mp.p)];
    n_p=numel(p)-1;
    logl=logl + log(p(1+ data.dx1));
else
    n_p=0;
end


if nargout >=3;  %% compute scores
    
    % step 1: compute derivative of contraction operator wrt. parameters
    dtdmp=zeros(X.n,1+ n_c + n_p);
    dtdmp(:, 1)=P*pk-1;
    dtdmp(:, 2:1+n_c)=-(P*dc).*pk; %
    if strcmp(method,'fullmle')
        vk=-cost+mp.beta*mp.ev;
        vr=-mp.RC-cost(1)+mp.beta*mp.ev(1);
        vmax=max(vk,vr);
        dtp=vmax+log(exp(vk-vmax)+exp(vr-vmax));

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

if nargout >=4;  %% compute Hessian

    %Hessian is computed FOR STRUCTURAL PARAMETERS ONLY
    if strcmp(method,'fullmle') || n_c>1
        error 'RTFM!'
    end

    % Parameters are: {RC, c}
    % Derivative of the diff between V(replace) and V(keep)
    % d(diff)/dRC
    ddiff(:,1)=-1+mp.beta*devdmp(1,1)-mp.beta*devdmp(:,1);
    % d(diff)/dc
    ddiff(:,2)=-dc(1)+mp.beta*devdmp(1,2)+dc-mp.beta*devdmp(:,2);

    % Derivative of pk
    % d(diff)/dRC
    dpk(:,1)=pk.*(pk-1).*ddiff(:,1);
    % d(diff)/dc
    dpk(:,2)=pk.*(pk-1).*ddiff(:,2);

    % Second derivatives of contraction mapping wrt parameters, RC and c
    d2t_dpi_dpj=nan(X.n,2,2);
    % (1,1) RC RC
    d2t_dpi_dpj(:,1,1) = P*dpk(:,1);
    % (1,2) RC c
    d2t_dpi_dpj(:,1,2) = P*dpk(:,2);
    % (2,1) c RC (assuming second deriv of c =0)
    d2t_dpi_dpj(:,2,1) = -P*( dc.*dpk(:,1) + dc(1)*(1-dpk(:,1)) );
    % (2,2) c c (assuming second deriv of c =0)
    d2t_dpi_dpj(:,2,2) = -P*( dc.*dpk(:,2) + dc(1)*(1-dpk(:,2)) );


    % Cross derivative of contraction mapping wrt EV and parameters RC and c
    d2tev_dev_dpi=nan(X.n,X.n,2);
    % wrt RC
    tmp=P(:,2:X.n).*repmat(dpk(2:X.n,1)',X.n,1);
    d2tev_dev_dpi(:,:,1)=mp.beta*[-sum(tmp,2) tmp];
    % wrt c
    tmp=P(:,2:X.n).*repmat(dpk(2:X.n,2)',X.n,1);
    d2tev_dev_dpi(:,:,2)=mp.beta*[-sum(tmp,2) tmp];

    % Second derivative of contracrtion fixed point, d2ev_dpi_dpj
    % Step 1: Computed components
    q=nan(X.n,2,2);
    for i=1:2
        for j=1:2
            q(:,i,j)=d2t_dpi_dpj(:,i,j)+d2tev_dev_dpi(:,:,i)*devdmp(:,j)+ ...
                                        d2tev_dev_dpi(:,:,j)*devdmp(:,i);
        end
    end

    % Step 2: put together (separately because have to cross-reference)
    for i=1:2 %by parameter
        d2ev_dpi_dpj(:,:,i) = F \ q(:,:,i);
    end

    %Compute Hessian of likelihood [2x2]
    % H=nan(2,2);
    % for i=1:2
    %     for j=1:2
    %         H(i,j)=sum( ...
    %             pk(data.x).*(pk(data.x)-1).*ddiff(data.x,i).*ddiff(data.x,j) + ...
    %             mp.beta*(lp-1+(data.d==2)).*(d2ev_dpi_dpj(1,i,j)-d2ev_dpi_dpj(data.x,i,j)) ...
    %                   );
    %     end
    % end

    H=nan(2,2);
    pkpr=pk(data.x).*(pk(data.x)-1);
    H(1,1)=sum(pkpr.*ddiff(data.x,1).*ddiff(data.x,1) + mp.beta*(lp-1+(data.d==2)).*(d2ev_dpi_dpj(1,1,1)-d2ev_dpi_dpj(data.x,1,1)));
    H(2,2)=sum(pkpr.*ddiff(data.x,2).*ddiff(data.x,2) + mp.beta*(lp-1+(data.d==2)).*(d2ev_dpi_dpj(1,2,2)-d2ev_dpi_dpj(data.x,2,2)));
    % H(1,2)=sum(pkpr.*ddiff(data.x,1).*ddiff(data.x,2) + mp.beta*(lp-1+(data.d==2)).*(d2ev_dpi_dpj(1,1,2)-d2ev_dpi_dpj(data.x,1,2)));
    H(1,2)=sum(pkpr.*ddiff(data.x,2).*ddiff(data.x,1) + mp.beta*(lp-1+(data.d==2)).*(d2ev_dpi_dpj(1,2,1)-d2ev_dpi_dpj(data.x,2,1)));
    H(2,1)=H(1,2);
    end


end % end of ll