% NFXP: procedure to compute contraction operator
%
%       NFXP.bellman :
%           Procedure to compute bellman equation
%
% See also: nfxp.bellman, nfxp.dbellman, nfxp.solve, nfxp.maxlik

classdef nfxp
    methods (Static)
        function [data] = simdata(X, mp, S, P, pk)
            N=S.N;
            T=S.T;
            
            if nargin ==3  % solve model if solution is not provided
                fxp_opt.printfxp=0;
                P = nfxp.statetransition(mp.p, X.n);
                cost=0.001*mp.c*X.grid;
                [ev, pk]=nfxp.solve(0, P, cost, mp, fxp_opt);
            end
            
            id=repmat((1:N),T,1);
            t=repmat((1:T)',1,N);
            u_dx=rand(T,N);
            u_d=rand(T,N);
            
            csum_p=cumsum(mp.p);
            dx1=0;
            for i=1:numel(csum_p);
                dx1=dx1+ (u_dx>(csum_p(i)));
            end;
            
            x=nan(T, N);
            x1=nan(T, N);
            x(1,:)=ones(1, N); % Intial conditions
            
            for it=1:T;
                d(it,:)=(u_d(it,:)<(1-pk(x(it,:)')'))+1;  % Keep =1, Replace =2
                x1(it,:)=min(x(it,:).*(d(it,:)==1) +(d(it,:)==2) + dx1(it,:), X.n);
                if it<T;
                    x(it+1,:)=x1(it,:);
                end
            end
            
            data.id=id;
            data.t=t;
            data.d=d;
            data.x=x;
            data.x1=x1;
            data.dx1=dx1;
            
            pfields=fieldnames(data);
            for i=1:numel(pfields);
                data.(pfields{i})=reshape(data.(pfields{i}), T*N,1);
            end
            
            
            
        end
        
        function data = readdata(X, bustypes)
            %Read data
            dta=load('busdata.txt','-ascii');
            disp('got here');
            % Cols in busdata1234:
            %     1: busid,
            %     2: bustype
            %     3: year
            %     4: month
            %     5: replacement_dummy
            %     6: accumulated mileage, x
            %     7: lagged accumulated milage, x1
            
            % Select busses
            if nargin>1
                dta=dta(dta(:,2) <=bustypes,:);
            end
            
            % Populate data struct with selected busdata
            data.id=dta(:,1);
            data.d=dta(:,5)+1; % d=1: keep engine; d=2: replace engine
            
            data.x=dta(:,6)/1000;     % Discretized accumulated mileage, t
            data.x1=dta(:,7)/1000;    % Discretized accumulated mileage, t-1
            data.dx1=data.x1-(data.d==1).*data.x;   % Discretized monthly mileage
        end
        
        function [data, grid] = discretize(data, X)
            % Discrete dx and compute x from discretized x
            discretize = @(x) ceil(x*X.n/X.max);
            data.x=discretize(data.x);       % Discretized accumulated mileage, t
            data.x1=discretize(data.x1);     % Discretized accumulated mileage, t-1
            data.dx1=(data.x1-(data.d==1).*data.x);   % Discretized monthly mileage
        end
        
        function P = statetransition(p, n)
            p=[p; (1-sum(p))];
            P=0;
            for i=0:numel(p)-1;
                P=P+sparse(1:n-i,1+i:n,ones(1,n-i)*p(i+1), n,n);
                P(n-i,n)=1-sum(p(1:i));
            end
            P=sparse(P);
        end
        
        
        function [ev]=E(v, p)
            ev=p(1)*v;
            for i=1:numel(p)-1
                ev=ev+ p(i+1)*[v(1+i:end); repmat(v(end), i, 1)];
                
            end
        end
        
        
        function [ev, pk]=bellman(ev, P, c, mp)
            % NFXP.BELMANN:     Procedure to compute bellman equation
            %
            % Syntax :          ev=nfpx.bellman(ev, P, c, mp);
            %
            % Inputs
            %  ev0              n x 1 matrix of expected values given initial
            %                   guess on value function
            % Output:
            %  ev1              n x 1 matrix of expected values given initial
            %                   guess of ev
            % See also:
            %   nfxp.maxlik, nfxp.solve, nfxp.dbellman
            V=-c+mp.beta*ev- (-mp.RC-c(1)+mp.beta*ev(1));
            maxV=max(V, 0);
            ev=P*(maxV + log(exp(V-maxV)  +  exp(-maxV)));
%           ev=nfxp.E(maxV + log(exp(V-maxV)  +  exp(-maxV)),mp.p);
            if nargout>1
                pk=1./(1+exp(-V));
            end
            
        end % end of NFXP.bellman
                
        
        function dev=dbellman(ev, pk,  P, mp)
            % NFXP.BELMANN:     Procedure to compute bellman equation
            %
            % Syntax :          ev=nfpx.dbellman(ev, pk,  P, mp);
            %
            % See also:
            %   nfxp.maxlik, nfxp.solve, nfxp.bellman
            
            n=numel(ev);
            
            tmp=P(:,2:n).*repmat(pk(2:n,1)',n,1);
            dev=(mp.beta*[-(sum(tmp,2)) tmp]);
            
            %Faster for large scale problems
            %tmp=bsxfun(@times, P(:,2:n), pk(2:n,1)');
            %?dev=sparse(mp.beta*[-(sum(tmp,2)) tmp]);
            
        end % end of NFXP.bellman
        
        
        function [ev, pk, F]=solve(ev, P, cost, mp, options)
            %NFXP.SOLVE: Contraction mapping fixed point poly-algorithm.
            %
            %  syntax:	[ Pk,EV,dEV]=nfxp.solve(ev, P, cost, mp, options):
            %
            %  INPUT:
            %     EV :      m x nc matrix or scalar zero. Initial guess of choice specific expected value function, EV.
            %     P:        State transition matrix
            %     cost:     State transition matrix
            %     mp:       Structure containing model parameters: mp.beta, mp.RC, mp.c, mp.ev
            %     options:  Optional setting for Fixed Point Algorithm
            %
            %     options.max_fxpiter (integer);
            %       Maximum number of times to switch between Newton-Kantorovich iterations and contraction iterations (default: 2).
            %
            %     options.min_cstp (integer):
            %       Minimum number of contraction steps (default: 6)
            %
            %     options.max_cstp (integer):
            %       Maximum number of contraction steps (default: 30)
            %
            %     options.ctol:
            %       Tolerance before switching to N-K algorithm (default: 0.1)
            %
            %     options.rtol:
            %       Relative tolerance before switching to N-K algorithm (default: 0.1)
            %
            %     options.nstep
            %       Maximum number of Newton-Kantorovich steps
            %
            %     options.ltol0
            %       Final exit tolerance in fixed point algorithm, measured in units of numerical precision (default is 1.0e-12)
            %
            %     options.printfxp
            %       Print iteration info for fixed point algorithm if > 0 (default is 0)
            %
            %  OUTPUT:
            %     Pk:     m x 1 matrix of conditional choice probabilities, P(d=keep|x)
            %     EV:     m x 1 matrix of expected value functions, EV(x)
            %     F:      m x m matrix of derivatives of Identity matrix
            %             minus contraction mapping operator, I-T' where T' refers to derivative of the expected  value function
            %
            % See also:
            %   nfxp.maxlik, nfxp.bellman
            
            %-----------------------------------------------------------------------------------------------------------------------------
            % Default settings for Fixed point algorithm
            
            global BellmanIter NKIter;

            opt.max_fxpiter= 1;             % Maximum number of times to switch between Newton-Kantorovich iterations and contraction iterations.
            opt.min_cstp   = 2;             % Minimum number of contraction steps
            opt.max_cstp   = 20;            % Maximum number of contraction steps
            opt.ctol       = 0.02;          % Tolerance before switching to N-K algorithm
            opt.rtol       = 0.02;          % Relative tolerance before switching to N-K algorithm
            opt.nstep      = 20;            % Maximum number of Newton-Kantorovich steps
            opt.ltol0      = 1.0e-12;       % Final exit tolerance in fixed point algorithm, measured in units of numerical precision
            opt.printfxp   = 1;             % print iteration info for fixed point algorithm if > 0
            
            if nargin>4
                pfields=fieldnames(options);
                for i=1:numel(pfields);
                    opt.(pfields{i})=options.(pfields{i});
                end
            end
            
            n=numel(cost);
            
            NKIter=0;
            for k=1:opt.max_fxpiter;
                % SECTION A: CONTRACTION ITERATIONS
                if opt.printfxp>0
                    fprintf('\n');

                    fprintf('Begin contraction iterations\n');
                    fprintf('   j           tol        tol(j)/tol(j-1) \n');
                end;
                
                tolm1=1;
                tic;
                
                for j=2:2:opt.max_cstp;
                    ev1=nfxp.bellman(ev, P, cost, mp);
                    ev=nfxp.bellman(ev1, P, cost, mp);
                    
                    BellmanIter=j;
                    
                    tolm=max(max(abs(ev-ev1)));
                    rtolm=tolm/tolm1;
                    
                    tolm1=tolm;
                    if opt.printfxp>0
                        fprintf(' %3.0f   %16.8f %16.8f\n',j, tolm,rtolm);
                    end
                    if tolm < opt.ctol;
                        break;
                    end;
                    if (j>=opt.min_cstp) && (abs(mp.beta*mp.beta-rtolm) < opt.rtol);
                        break
                    end
                end
                
                % SECTION 2: NEWTON-KANTOROVICH ITERATIONS
                if opt.printfxp>0
                    fprintf('\n');
                    fprintf('Begin Newton-Kantorovich iterations at time: %1.5f\n', toc);
                    fprintf('  nwt          tol   \n');
                    
                end
                [ev1, pk]=nfxp.bellman(ev, P, cost, mp);
                for nwt=1:opt.nstep;
                    F=speye(n) - nfxp.dbellman(ev, pk, P,  mp) ;
                    
                    ev1=ev-F\(ev-ev1);
                    
                    [ev, pk]=nfxp.bellman(ev1, P, cost, mp);
                    tolm=max(max(abs(ev-ev1)));
                    
                    if opt.printfxp>0
                        fprintf('   %d     %16.8e  \n',nwt,tolm);
                    end
                    adj=ceil(log(abs(max(max(ev)))));
                    ltol=opt.ltol0*10^adj;  % Adjust final tolerance
                    if (tolm < ltol) || (nwt > opt.nstep);
                        break
                    end
                    temp=ev1; % FIXME - Whats going on here
                    ev1=ev;
                    ev=temp;
                    NKIter=nwt;
                end; % FXP iterations --- holds both the contraction and N-K iterations
                
                if (tolm < ltol);
                    break; % convergence achieved in N-K iterations
                end;
                
                if opt.printfxp>0 && NKIter < opt.nstep;
                    fprintf('\n');
                    fprintf('Time to convergence is %10.8f seconds\n', toc);
                elseif NKIter >= opt.nstep;
                    warning('Maximum number of iterations exceeded without convergence!');
                end;
                
            end % end of fxp iterations.
        end % end of solve
        
        
        function [mphat, mpse, cov, g, ll, iterinfo]=maxlik(llfun, p0, pnames, options);
            %NFXP.MAXLIK:  Routine to maximize likelihood functions based on
            %               - Analytical gradients.
            %               - BHHH Algorithm
            %               - Line search (step-halfing)
            %
            %   Syntax:   [mphat, mpse, cov, g, ll, convflag]=maxlik(llfun, p0, pnames, options);
            %
            %
            %   OUTPUTS:
            %       mphat:      Estimated parameters (structure with as many fields as p0)
            %       mpse:       Standard errors of estimated parameters (structure with number of fields equal to number elements in pnames, nP)
            %       cov:        Variance Covariance Matrix (nP x nP matrix)
            %       g:          Score vectors (N x nP)
            %       ll:         Likelihood contributions (N x 1 vector)
            %       convflag    Bolean (true/false) idicator of convergence
            %
            %   INPUTS:
            %       llfun:      Function handle for likelihood function
            %
            %       p0:         Structure holding starting values of parameters
            %
            %       pnames:     k dimensional cell array, that holds names of the field in the parameter structure
            %                   to be estimated. Then likelihood is only maximized with respect to the parameters in pnames
            %                   and it is thus possible to fix parameters at the values set in p0
            %
            %                   Example:
            %                       Suppose parameter structure contains both p0.delta0 and p0.delta1
            %                       If pnames={'delta0','delta1'}, then both delta0 and delta1 are estimated
            %                       If pnames={'delta0'}, then only delta0 is estimated and delta1 is fixed at starting value (set in p0)
            %
            %
            %       options:    Structure holding voluntary options for ML. If unspecified default options are used
            %                   options structure has the following fields:
            %
            %           options.title (String):
            %                   A title for the model.
            %
            %           options.output (Scalar) :
            %                   0 if no output, 1 for final output, 2 for detailed iteration output (default is 2)
            %
            %           options.maxit (Scalar):
            %                   Maximum number of iterations (default is 100)
            %
            %           options.tolg (Scalar):
            %                   Tolerance for gradient'direc (default is 1e-3)
            %
            %           options.lambda0:
            %                   The basic step length (default is 1).
            % See also:
            %   nfxp.solve, nfxp.bellman
            
            tt=tic;
            par.title={''};
            par.maxit=100;
            par.tolg=1e-1;
            par.output=1;
            par.lambda0=1;
            par.method='Partial MLE';
            p1=p0;

            global BellmanIter NKIter
            iterinfo.MajorIter=0;
            iterinfo.ll=0;
            iterinfo.NKIter=0;
            iterinfo.BellmanIter=0;


                        
            if nargin==4
                if isfield(options, 'title');
                    par.title=options.title;
                end
                if isfield(options, 'output');
                    par.output=options.output;
                end
                if isfield(options, 'maxit');
                    par.maxit=options.maxit;
                end
                if isfield(options, 'tolg');
                    par.tolg=options.tolg;
                end
                if isfield(options, 'lambda0');
                    par.lambda0 = options.lambda0;
                end;
                if isfield(options, 'method');
                    par.methods = options.method;
                end;
            end
            
            lambda=par.lambda0;   % initialize stepsize
            k = length(pnames);  % Number of parameters to estimate
            
            % -------------------------
            %  ** BHHH Optimization **
            % -------------------------

            [ll0,ev,s0]=llfun(p0);
            iterinfo.ll=iterinfo.ll+1;
            iterinfo.NKIter=iterinfo.NKIter+NKIter;
            iterinfo.BellmanIter=iterinfo.BellmanIter+BellmanIter;


            h=s0'*s0;
            g=sum(s0)';
            direc=h\g;
            tolm=g'*direc;

            iterinfo.MajorIter=0;
            for it=0:par.maxit;
                for i=1:k;
                    np=numel(p1.(char(pnames(i))));
                    p1.(char(pnames(i)))=p0.(char(pnames(i)))+ lambda*direc(i:i+np-1);
                    
                end
                p1.ev=ev;
                
             
                if lambda==par.lambda0;   % If previous step was accepted
                    [ll,ev,s]=llfun(p1);
                    iterinfo.ll=iterinfo.ll+1;
                    iterinfo.NKIter=iterinfo.NKIter+NKIter;
                    iterinfo.BellmanIter=iterinfo.BellmanIter+BellmanIter;


                else
                    [ll, ev, s]=llfun(p1); % Do we need to recalculate s, when line searching
                    iterinfo.ll=iterinfo.ll+1;
                    iterinfo.NKIter=iterinfo.NKIter+NKIter;
                    iterinfo.BellmanIter=iterinfo.BellmanIter+BellmanIter;


                end
                if sum(ll)<sum(ll0);
                    lambda=lambda/2;
                    if par.output>1;
                        fprintf('\n');
                        fprintf('.............LINESEARCHING\n');
                        fprintf('\n');
                    end
                    if lambda<=.01;
                        if par.output>1;
                            fprintf('.............WARNING: Failed to increase - downhill step accepted\n');
                        end
                        p0=p1;
                        ll0=ll;                 % Accept step
                        lambda=par.lambda0;
                    end
                else       % IF INCREASE IN LIKELIHOOD, CALCULATE NEW DIRECTION VECTOR
                    p0=p1;
                    ll0=ll;                 % Accept step
                    
                    
                    % ** Plot iteration info ***
                    if par.output>1;
                        fprintf('Iteration: %d   \n', iterinfo.MajorIter);
                        fprintf('grad*direc       %10.5f \n', g'*direc);
                        fprintf('Log-likelihood   %10.5f \n', sum(ll));
                        fprintf('Step length      %10.5f \n', lambda);
                        fprintf('---------------------------------------------------------------------\n');
                        fprintf('    Param.       Estimates      Direction         Gradient\n');
                        i=0;
                        for iP=1:numel(pnames);
                            for j=1:numel(p0.(char(pnames(iP))));
                                i=i+1;
                                fprintf('    %-10s %10.4f    %10.4f        %10.4f\n', char(pnames(iP)), p0.(char(pnames(iP)))(j), direc(i), g(i));
                            end;
                        end;
                        fprintf('---------------------------------------------------------------------\n');

                        fprintf('\n');
                    end
                    h=s'*s;
                    g=sum(s)';
                    iterinfo.MajorIter=iterinfo.MajorIter+1;
                    direc=h\g;
                    tolm=g'*direc;
                    lambda=par.lambda0; % and reset stepsize
                end
               
                if g'*h*g < 1;  % Stopping rule
%                    g'*h*g
                    break;
                end
 
                
                if tolm < par.tolg;  % Stopping rule
                    break;
                end
                
            end
            
            cov=inv(h); 
            se=(sqrt(diag(cov)));
            mphat=p0;
            i=0;
            for iP=1:numel(pnames);
                for j=1:numel(p0.(char(pnames(iP))))
                    i=i+1;
                    mpse.(char(pnames(iP)))(j)=se(i);
                end
            end
            
            % ---------------------------
            %  ** Plot final output  ***
            % ---------------------------
            if it<par.maxit;
                if par.output >= 1;
                    cr=corrcov(cov);
                    fprintf('---------------------------------------------------------------------\n');
                    fprintf('***                   Convergence Achieved                        ***\n');
                    fprintf('---------------------------------------------------------------------\n');
                    disp('                         _ ');
                    disp('                         \`\ ');
                    disp('                         |= | ');
                    disp('                        /-  ;.---. ');
                    disp('                  _ __.''     (____) ');
                    disp('                   `         (_____) ');
                    disp('                   _''  ._ .'' (____) ');
                    disp('                    `        (___) ');
                    disp('                   --`''------''` ');
                    fprintf('%s \n', char(par.title));
                    fprintf('Estimation Method: %s\n', par.method);
                    fprintf('Number of iterations: %d\n', iterinfo.MajorIter);
                    fprintf('grad*direc       %10.5f \n', g'*direc);
                    fprintf('Log-likelihood   %10.5f \n', sum(ll));
                    fprintf('Step length      %10.5f \n', lambda);
                    fprintf('\n');
                    fprintf('    %-13s %13s %13s %13s\n','Param.','Estimates','s.e.','t-stat');
                    fprintf('---------------------------------------------------------------------\n');
                     i=0;
                     for iP=1:numel(pnames);
                            for j=1:numel(p0.(char(pnames(iP))))
                                i=i+1;
                                myStr = char(pnames(iP));
                                if numel(p0.(char(pnames(iP))))>1; 
                                    myStr = sprintf('%-1s(%i)',char(pnames(iP)),j);
                                end;
                                fprintf('    %-13s %13.4f %13.4f %13.4f\n', myStr, p0.(char(pnames(iP)))(j), se(i), p0.(char(pnames(iP)))(j)/se(i));
                                mpse.(char(pnames(iP)))(j)=se(i);
                            end
                    end
                    fprintf('---------------------------------------------------------------------\n');

                    fprintf('\n');

                    fprintf('Correlation matrix of parameters\n');
                    fprintf('%18s','');
                    i=0;
                    for iP=1:k;
                        for j=1:numel(p0.(char(pnames(iP))))
                            i=i+1;
                            myStr = char(pnames(iP));
                            if numel(p0.(char(pnames(iP))))>1; 
                                myStr = sprintf('%-s(%i)',char(pnames(iP)),j);
                            end;

                            fprintf('%18s', myStr);
                        end;
                    end
                    fprintf('\n');

                    i=0;
                    for iP=1:k;
                        for j=1:numel(p0.(char(pnames(iP))))
                            i=i+1;
                            myStr = char(pnames(iP));
                            if numel(p0.(char(pnames(iP))))>1; 
                                myStr = sprintf('%-1s(%i)',char(pnames(iP)),j);
                            end;

                            fprintf('%-18s', myStr);
                            for j2=1:size(cr,1);
                                fprintf('%18.4f',cr(i,j2));
                            end
                            fprintf('\n');
                        end;
                    end
                    if ismac & par.output >= 3 
                        !say 'Success! Convergence Achieved'
                    end
                    tt_end=toc(tt);
                    fprintf('\n');
                    fprintf('Time to convergence is %3.0f min and %2.2f seconds\n', floor((tt_end)/60),tt_end - 60*floor((tt_end)/60));

                    fprintf('\n');
                    
                end
                % output a flag that we did converge
                iterinfo.convflag = true;
            else
                fprintf('----------------------------------------------------------------------------\n');
                fprintf('***   BHHH failed to converge: Maximum number of iterations is reached *** \n');
                fprintf('----------------------------------------------------------------------------\n');
                disp('                   _,....._ ')
                disp('                  (___     `''-.__ ')
                disp('                 (____ ')
                disp('                 (____ ')
                disp('                 (____         ___ ')
                disp('                      `)   .-''` ')
                disp('                      /  .'' ')
                disp('                     | =| ')
                disp('                      \_\ ')
                
                if ismac
                    !say 'B triple H failed to converge: Maximum number of iterations is reached without convergence. Better luck next time'
                end
                
                % no luck this time
                iterinfo.convflag = false;
                
            end

        end % end of estim.maxlik
        
    end % end of methods
    
    
end % end of estim class