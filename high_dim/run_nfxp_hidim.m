function out=run_nfxp_hidim(param,data,X0)	
%This funtion runs NFXP-NK estimator for the series of dataset
%
%INPUT: structure of common parameters
%				MC_* simulated data
%				X0 matrix of starting values (first 7 rows are used, columns=multistarts)
%OUTPUT: self-explaining structures of results for Rust87 implementation (golden standard) 
%				 and Judd-Su design replicated exactly

%Convert parameters to our notation
nMC=param.MC;
X.n=param.N;
X.grid= (0:X.n-1)';
S.N=param.nBus;
S.T=param.nT;

%parameters of algorithms
fxp_opt.printfxp=0;		% 0: No output,1: Fixed point iteration output 1   

ml_opt.title='Structural parameters in Rust(1987) Engine replacement model';
ml_opt.output=0;     % 0: No output, 1: Only final output,2: Detailed iteration output
ml_opt.method='Partial MLE';  % 'Full MLE' or 'Partial MLE'

out.Runs=0;
out.TotalSuccess=0;
out2.Runs=0;
out2.TotalSuccess=0;
out2.BellmanIter=0;
out2.NKIter=0;

%Monte Carlo loop
for j=1:param.MC

	%convert simulated data to Bertel's format
	datasim.d=reshape(data.MC_dt(:,:,j)+1,S.T*S.N,1);
	datasim.x=reshape(data.MC_xt(:,:,j),S.T*S.N,1);
	datasim.dx1=reshape(data.MC_dx(:,:,j)-1,S.T*S.N,1);

	%Only run golden standard (Rust 87)
	% Frequency estimator for transition probs
	mp_start.beta=param.beta;
	mp_start.RC=X0(2,j); %use JS initial values for the variables on the first stage
	mp_start.c=X0(1,j);
	mp_start.ev=0;
	mp_start.p=param.thetaProbs(1:end-1); %freq estimator
	P = nfxp.statetransition(mp_start.p, X.n);

  % Likelihood function wrappers
  ll_p=@(mp) ll(datasim, X, P, mp, fxp_opt, 'pmle');
	% Partial likelihood for RC and c (we are only doing this here)
	estimationtimer=tic;
  pnames={'RC','c'};
	[mphat, mpse, cov, g, llval, iterinfo1]=nfxp.maxlik(ll_p, mp_start, pnames, ml_opt);
	timetoestimate=toc(estimationtimer);
	%save results
	out.param=param;
	out.mc_iter(j)=j;
	out.Runs=out.Runs+1;
	out.TotalSuccess=out.TotalSuccess+iterinfo1.convflag;
	out.converged(j)=iterinfo1.convflag;
	%estimated parameters
	out.RC.est(j)=mphat.RC;
	out.RC.se(j)=mpse.RC;
	out.c.est(j)=mphat.c;
	out.c.se(j)=mpse.c;
	%performance indicators
	out.runtime(j)=timetoestimate;
	out.MajorIter(j)=iterinfo1.MajorIter;
	out.FuncEval(j)=iterinfo1.ll;
	out.BellmanIter(j)=iterinfo1.BellmanIter;
	out.NKIter(j)=iterinfo1.NKIter;
end

end %function



