clear 
close all
clf

notes=0;
load('Results_nBus');
tit='N=175'

N_ampl=nan(numel(ampl), 1);
cputime_ampl=nan(numel(ampl), 1);
N_nfxp=nan(numel(nfxp), 1);
cputime_nfxp=nan(numel(nfxp), 1);
for i=1:numel(ampl)
  N_ampl(i)=nBus(i)*120;
  cputime_ampl(i)=ampl{i}.run_time/ampl{i}.num_iter;
  N_nfxp(i)=nBus(i)*120;
  cputime_nfxp(i)=mean(nfxp{i}.runtime./nfxp{i}.MajorIter);
end

[B_nfxp,BINT,R,RINT,STATS_nfxp]=regress((cputime_nfxp),[ones(numel(N_nfxp),1) N_nfxp]);
[B_ampl,BINT,R,RINT,STATS_ampl]=regress((cputime_ampl),[ones(numel(N_ampl),1) N_ampl]);
N=120*[0; nBus'];

ms=4;
lw=2;
scrsz = get(groot,'ScreenSize');
figure('Name','CPU time per major iteration','Position',[1 scrsz(3)/3 scrsz(3)/3 scrsz(3)/3])
hold on
p_mpec=plot(N_ampl,([cputime_ampl]), 'ok', 'LineWidth',lw, 'MarkerSize',ms);
set(gca,'FontSize',14);
hold on
plot(N,B_ampl(1)+B_ampl(2)*(N),'--b', 'LineWidth',lw, 'MarkerSize',ms);


p_nfxp=plot(N_nfxp,([cputime_nfxp]), '*r', 'LineWidth',lw, 'MarkerSize',ms);
xlabel('Sample size');
ylabel('cpu time per major iteration (seconds)');

hold on
plot(N,B_nfxp(1)+B_nfxp(2)*(N), '-r', 'LineWidth',lw, 'MarkerSize',ms);
nfxp_fit= sprintf('%8.3f + %8.2f*NT/100.000 (R^2 = %8.3f)',B_nfxp(1), B_nfxp(2)*1e6, STATS_nfxp(1))
R2_nfxp=STATS_nfxp(1)

ampl_fit= sprintf('%8.3f + %8.2f*NT/100.000 (R^2 = %8.3f)',B_ampl(1), B_ampl(2)*1e6, STATS_ampl(1))
R2_ampl=STATS_ampl(1)


if notes==1;
title(tit);
legend('MPEC-AMPL', ampl_fit,  'NFXP-NK', nfxp_fit, 'location',  'NorthWest');
else
legend([p_mpec p_nfxp], {'MPEC-AMPL', 'NFXP-NK'}, 'location',  'NorthWest');
end
legend('boxoff');
xlim([0 nBus(end)*120])
axis 'square'



figure('Name','Total CPU time','Position',[1 scrsz(3)/3 scrsz(3)/3 scrsz(3)/3])
for i=1:numel(ampl)
  N_ampl(i)=nBus(i)*120;
  cputime_ampl(i)=ampl{i}.run_time;
  N_nfxp(i)=nBus(i)*120;
  cputime_nfxp(i)=mean(nfxp{i}.runtime);
end

[B_nfxp,BINT,R,RINT,STATS_nfxp]=regress((cputime_nfxp),[ones(numel(N_nfxp),1) N_nfxp]);
[B_ampl,BINT,R,RINT,STATS_ampl]=regress((cputime_ampl),[ones(numel(N_ampl),1) N_ampl]);

hold on
p_mpec=plot(N_ampl,([cputime_ampl]), 'ok', 'LineWidth',lw, 'MarkerSize',ms);
set(gca,'FontSize',14);

hold on
plot(N,B_ampl(1)+B_ampl(2)*(N),'--b', 'LineWidth',lw, 'MarkerSize',ms);


p_nfxp=plot(N_nfxp,([cputime_nfxp]), '*r', 'LineWidth',lw, 'MarkerSize',ms);
xlabel('Sample size');
ylabel('total cpu time (seconds)');

hold on
plot(N,B_nfxp(1)+B_nfxp(2)*(N), '-r', 'LineWidth',lw, 'MarkerSize',ms);
nfxp_fit= sprintf('E(total cpu time|Sample size) = %8.3f + %8.2f*NT/100.000 (R^2 = %8.3f)',B_nfxp(1), B_nfxp(2)*1e6, STATS_nfxp(1))
R2_nfxp=STATS_nfxp(1)

ampl_fit= sprintf('E(total cpu time|Sample size)= %8.3f + %8.2f*NT/100.000 (R^2 = %8.3f)',B_ampl(1), B_ampl(2)*1e6, STATS_ampl(1))
R2_ampl=STATS_ampl(1)


if notes==1;
title(tit);
legend('MPEC-AMPL', ampl_fit,  'NFXP-NK', nfxp_fit, 'location',  'NorthWest');
else
legend([p_mpec p_nfxp], {'MPEC-AMPL', 'NFXP-NK'}, 'location',  'NorthWest');
end
legend('boxoff');

xlim([0 nBus(end)*120])
axis 'square'


