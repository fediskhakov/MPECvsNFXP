function RTable(param,dimvec,ampl,nfxp,doplots)
	%this function makes nice tables of results from the MC run
	%it can be called during the MC loop to save intermediate results on disk

	%number of produced results for different betas
	numdims=min([size(dimvec,2);numel(ampl);numel(nfxp)]);
	
	fprintf('\n\n\n');
	diary (['GrandResultsMC_hidim' param.out_suff '.txt']);
	disp(datestr(now));

	fprintf('Table 1: Estimation results (RMSE in brackets)\n');
	fprintf('--------------------------------------------------------------\n');
	fprintf('%20s %18s %18s \n','Method','RC','c');
	fprintf('--------------------------------------------------------------\n');
	fprintf('%20s %18.4f %18.4f \n','True values -->     ',param.RC,param.thetaCost);
	fprintf('--------------------------------------------------------------\n');
	for i=1:numdims
	fprintf('N = %1.0f, M = %1.0f\n',dimvec(1,i),dimvec(2,i));
	fprintf('--------------------------------------------------------------\n');
	%MPEC/AMPL
	if ~param.dontAMPL
	fprintf('%20s %8.4f(%8.4f) %8.4f(%8.4f)\n',...
					'MPEC/AMPL',...
					ampl{i}.mean_estimates(2),ampl{i}.rmse_estimates(2),...
					ampl{i}.mean_estimates(1),ampl{i}.rmse_estimates(1)...
				 )
	end
	fprintf('%20s %8.4f(%8.4f) %8.4f(%8.4f)\n',...
					'NFXP-NK/BHHH',...
					mean(nfxp{i}.RC.est),sqrt(mean((nfxp{i}.RC.est-param.RC).^2)),...
					mean(nfxp{i}.c.est),sqrt(mean((nfxp{i}.c.est-param.thetaCost).^2))...
				 )
	fprintf('--------------------------------------------------------------\n');
	end
	fprintf('\n\n');

	fprintf('Table 2: Numerical performance\n');
	fprintf('-----------------------------------------------------------------------------------------------------------------------\n');

	fprintf('%14s  %14s %14s %14s %14s %14s %14s\n', ''    , 'Runs Converged', 'CPU Time' ,'# of Major'	,'# of Func.'	,'# of Contract','# of N-K');
	fprintf('%-14s %14s %14s %14s %14s %14s %14s\n', 'N, M', sprintf('(out of %g)',param.MC) ,'(in sec.)','Iter'			,'Eval.'		,'Iter.','Iter.') 
	if ~param.dontAMPL
	fprintf('-----------------------------------------------------------------------------------------------------------------------\n');
	fprintf('                                                  MPEC/AMPL\n');
	fprintf('-----------------------------------------------------------------------------------------------------------------------\n');
	for i=1:numdims
		fprintf('%6.0f, %-6.0f %14d %14.4f %14.4f %14.4f %14s %14s\n', ...
			dimvec(1,i),dimvec(2,i),  ampl{i}.TotalSuccess(1) ,mean(ampl{i}.runtime(ampl{i}.converged(:,1)==1,1)) , ...
			ampl{i}.num_iter(1),ampl{i}.ave_feval(1), '-','-') 
	end
	end %dontAMPL
	fprintf('-----------------------------------------------------------------------------------------------------------------------\n');
	fprintf('                                                  NFXP-NK with 2 step ML and BHHH\n');
	fprintf('-----------------------------------------------------------------------------------------------------------------------\n');
	for i=1:numdims
		fprintf('%6.0f, %-6.0f %14d %14.4f %14.4f %14.4f %14.4f %14.4f\n', ...
			dimvec(1,i),dimvec(2,i),  nfxp{i}.TotalSuccess ,mean(nfxp{i}.runtime(nfxp{i}.converged)) , ...
			mean(nfxp{i}.MajorIter),mean(nfxp{i}.FuncEval), ...
			mean(nfxp{i}.BellmanIter),mean(nfxp{i}.NKIter))
	end
	fprintf('-----------------------------------------------------------------------------------------------------------------------\n');
	fprintf('\n');
	fprintf('Remaining parameters\n');
	fprintf('RC   = %10.5f \n',param.RC);
	fprintf('c0   = %10.5f \n',param.thetaCost);
	for i= 1:numel(param.thetaProbs)
		fprintf('p(%d) = %10.5f \n',i, param.thetaProbs(i));
	end
	fprintf('n    = %10.5f \n',param.N);
	fprintf('N    = %10.5f \n',param.nBus);
	fprintf('T    = %10.5f \n',param.nT);
	fprintf('beta = %10.5f \n',param.beta);
	fprintf('True parameters of milage transition: \n');
	fprintf('alfa = %10.5f \n',param.jr.alpha);
	fprintf('sigma= %10.5f \n',param.jr.sigma_dx);
	fprintf('bound= %10.5f \n',param.jr.upperbnd);
	fprintf('\n\n');

	diary off

	if doplots && ~param.dontAMPL
		for i=1:numdims
			fig=figure('Color',[1 1 1],'NextPlot','new');
			ax=axes('Parent',fig);
			xx=(0.001:0.001:1)';
			data{1}=ampl{i}.runtime(ampl{i}.converged(:,1)==1,1);
			hold on
			plot(quantile(data{1},xx),xx, '-k','DisplayName','MPEC/AMPL','Parent',ax);
			% plot(quantile(ampl{i}.runtime(:,2),xx),xx, ':k','DisplayName','MPEC/ktrlink','Parent',ax);
			% plot(quantile(ampl{i}.runtime(:,3),xx),xx, ':k','DisplayName','NFXP-SA with ktrlink','Parent',ax);
			% plot(quantile(result_nk{i}.runtime,xx),xx, ':k','DisplayName','NFXP-NK with ktrlink','Parent',ax);
			plot(quantile(nfxp{i}.runtime,xx),xx, '-r','DisplayName','NFXP-NK with 2 stage ML and BHHH','Parent',ax);
			xlabel('CPU time (seconds)')
			legend1 = legend(ax,'show');
			set(legend1,'EdgeColor',[1 1 1],'Location','SouthEast','YColor',[1 1 1],'XColor',[1 1 1]);
			title(sprintf('Distribution of estimation CPU Time with N=%1.0f, M=%1.0f (conditional on converging)',dimvec(1,i),dimvec(2,i)));
		end
	end
end