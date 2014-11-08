function [c,ceq,DC,DCeq] = confuneq(X)

% This m-file computes the value of the constraints (c, ceq) 
% and its Jacobian (DC, DCeq) in the constrained optimization formulation 
% of the Harold Zurcher bus-engine replacement model.
%
% source: Su and Judd (2011), Constrained Optimization Approaches to
% Estimation of Structural Models.
% Code Revised: Che-Lin Su, May 2010.
% Code Revised: Iskhakov, Rust, Schjerning, October 2014

global N M beta CbEVPrime TransProbPrime RPrime EVPrime;

[CbEV, EV, TransProb, thetaProbs, RC]  = DefVar(X);

%
% Bellman equation for states below N-M
indices = repmat((1:N)',1,M)+repmat((1:M),N,1)-1;

CbEV(N+1:N+M)=CbEV(N);
v1=CbEV(indices);
v2=-RC+CbEV(1);
vmax=max(v1,v2);
s1 = exp(v1-vmax);
s2 = exp(v2-vmax);
s =  s1 + s2;
logs = vmax+log(s);
BellmanViola = logs * TransProb - EV;

% Defne and evaluate nonlinear inequality constraints
c = [];
% Define and evaluate nonlinear equality constraints
ceq = BellmanViola;

% Define and evaluate the constraint Jacobian (DC, DCeq).   
    if nargout > 3
        DC= [];
    
        %d1 = ((CbEVPrime.*repmat(exp(CbEV),1,length(X)) + exp(-RC+CbEV(1))*repmat(RPrime,N+M,1)))./(repmat(exp(CbEV)+exp(-RC+CbEV(1)),1,length(X)));
        maxCbEV=max([CbEV;-RC+CbEV(1)]);%scalar
        d1 = ((CbEVPrime.*repmat(exp(CbEV-maxCbEV),1,length(X)) + exp(-RC+CbEV(1)-maxCbEV)*repmat(RPrime,N+M,1)))./(repmat(exp(CbEV-maxCbEV)+exp(-RC+CbEV(1)-maxCbEV),1,length(X)));
        sum1 =  reshape(sum(reshape(d1(indices',:) .* repmat(TransProb,N,length(X)),M,N, length(X) )),N,length(X)); 
        sum2 = logs*TransProbPrime';
        DCeq = sum1 + sum2 - EVPrime;
        DCeq = DCeq' ;
    end
end
