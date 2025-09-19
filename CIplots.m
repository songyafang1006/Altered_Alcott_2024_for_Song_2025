
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Plotting script   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Need to transpose initially to get Confidence intervals
d13CCI = real(sens.Aiso.');
CO2CI = real(sens.CO2atm.');
O2ACI = real(sens.O2_A.');
O2DPCI = real(sens.O2_DP.');
Dist_Preac_BurialCI = real(sens.Dist_Preac_Burial.');
SRP_DPCI = real(sens.SRP_DP.');
fanoxicdistCI = real(sens.fanoxicdist.');
GASTCI = real(sens.GAST.');
forgCI = real(sens.forg.');

%% Confidence intervals
d13C95quant = quantile(d13CCI,[0.05 0.95]);
CO295quant = quantile(CO2CI,[0.05 0.95]);
O2A95quant = quantile(O2ACI,[0.05 0.95]);
O2DP95quant = quantile(O2DPCI,[0.05 0.95]);
Dist_Preac_Burial95quant = quantile(Dist_Preac_BurialCI,[0.05 0.95]);
SRP_DP95quant = quantile(SRP_DPCI,[0.05 0.95]);
fanoxicdist95quant = quantile(fanoxicdistCI,[0.05 0.95]);
GAST95quant = quantile(GASTCI,[0.05 0.95]);
forg95quant = quantile(forgCI,[0.05 0.95]);

%% Medians
d13C95median = quantile(d13CCI,[0.5]);
CO295median = quantile(CO2CI,[0.5]);
O2A95median = quantile(O2ACI,[0.5]);
O2DP95median = quantile(O2DPCI,[0.5]);
Dist_Preac_Burial95median = quantile(Dist_Preac_BurialCI,[0.5]);
SRP_DP95median = quantile(SRP_DPCI,[0.5]);
fanoxicdist95median = quantile(fanoxicdistCI,[0.5]);
GAST95median = quantile(GASTCI,[0.5]);
forg95median = quantile(forgCI,[0.5]);

global starting
%%%%%% define colours
c_mean = [0.2 0.6 0.6] ;
c_std = [0.3 0.7 0.7] ;
c_range = [ 0.4 0.8 0.8] ;

%%%% output to screen
fprintf('running sens plotting script... \t')
tic

%%%% make column vector
sens.time_myr = sens.time(:,1) /1e6 ;
%% Call in data
figure('Color',[0.80 0.80 0.70])

%d13C
load('Havigd13C.mat')
d13Ctime = -1 * Havigd13C(:,1);
d13Ccalcite = Havigd13C(:,2);
d13Cdolomite = Havigd13C(:,3);
d13Cother = Havigd13C(:,4);

subplot(4,2,1)
scatter(d13Ctime, d13Ccalcite,5,'o')
hold on
scatter(d13Ctime, d13Cdolomite,5,'+')
hold on
scatter(d13Ctime, d13Cother,5,'*')
hold on
plot((sens.time_myr),(d13C95median),'linewidth',3,'color',c_mean)
hold on
plot((sens.time_myr),(d13C95quant),'linewidth',1,'color',c_range)

xlim([-4e3 0]);
ylim([-25 20])
legend('calcite','dolomite','other')

%CO2
load('CO2proxy')


subplot(4,2,2)
box on
xlabel('Time (Ma)')
ylabel('CO2')
%%%% plot this model
semilogy((sens.time_myr),(CO295median*1e6),'linewidth',3,'color',c_mean)
hold on
semilogy((sens.time_myr),(CO295quant*1e6),'linewidth',1,'color',c_range)
hold on
semilogy(1e3*HighTime, HighCO2*1e6,'linewidth',2,'color','k')
hold on
semilogy(1e3*LowTime, LowCO2*1e6,'linewidth',2,'color','k')
xlim([-4e3 0]);
ylim([10 1e7]);

title('CO2')

%Atmos O2
load('AtmosO2proxy.mat')
O2_A_min = min((sens.O2_A/3.7e19),[],2) ;
O2_A_max = max((sens.O2_A/3.7e19),[],2) ;

subplot(4,2,3)
box on
xlabel('Time (Ma)')
ylabel('PAL O2')
%%%% plot this model
semilogy((sens.time_myr),(O2A95quant/3.7e19),'linewidth',1,'color',c_range)
hold on
semilogy((sens.time_myr),(O2A95median/3.7e19),'linewidth',3,'color',c_mean)
hold on
semilogy(AtmosO2proxtime, AtmosO2proxlow,'linewidth',2,'color','k')
hold on
semilogy(AtmosO2proxtime, AtmosO2proxhigh,'linewidth',2,'color','k')
xlim([-4e3 0]);
ylim([1e-7 10])
title('O2 atmosphere PAL')

%O2 Deep
subplot(4,2,4)
box on
semilogy((sens.time_myr),(O2DP95median/2.21e17),'linewidth',3,'color',c_mean)
hold on
semilogy((sens.time_myr),(O2DP95quant/2.21e17),'linewidth',1,'color',c_range)
xlim([-4e3 0]);
ylim([1e-9 10])
xlabel('Time (Ma)')
ylabel('deep O2 relative')
title('Relative O2 Deep')

%P burial
load('ReinhardPdata.mat')

XXX = [Time Reinhard_P] ;
ZZ = XXX(~any(isnan(XXX)| isinf( XXX ), 2 ),: ) ;
x = ZZ(:,1) ;
y = ZZ(:,2) ;
edges = (0:50:4e3);
[~,~,loc]=histcounts(x,edges);
meany = accumarray(loc(:),y(:))./accumarray(loc(:),1);
xmid = 0.5*(edges(1:end-1)+edges(2:end));
xmid = xmid.' ;
xmid = xmid(1:end-10,:);

subplot(4,2,5)
box on
xlabel('Time (Ma)')
yyaxis right
scatter((-1*Time), logReinhard_P,5)
hold on
plot((-1*xmid),log10(meany),'linewidth',2,'color','k')
ylim([-3 1])

yyaxis left
set(gca, 'YScale', 'log')
ylabel('P burial')
plot((sens.time_myr),(Dist_Preac_Burial95median),'linewidth',3,'color',c_mean)
hold on
plot((sens.time_myr),(Dist_Preac_Burial95quant),'linewidth',1,'color',c_range)
ylim([1e7 1e13])

xlim([-4e3 0]);
title('Total P burial')

%P deep
subplot(4,2,6)
box on
semilogy((sens.time_myr),(SRP_DP95median/2790e12),'linewidth',3,'color',c_mean)
hold on
semilogy((sens.time_myr),(SRP_DP95quant/2790e12),'linewidth',1,'color',c_range)
xlim([-4e3 0]);
ylim([1e-3 5])
title('P deep relative')


%fanoxic
subplot(4,2,7)
plot((sens.time_myr),(fanoxicdist95median),'linewidth',3,'color',c_mean)

hold on
plot((sens.time_myr),(fanoxicdist95quant),'linewidth',1,'color',c_range)
title('fanoxicdist')
xlim([-4e3 0]);
ylim([0 1]);

%Temperature
subplot(4,2,8)
plot((sens.time_myr),(GAST95median-273),'linewidth',3,'color',c_mean)
hold on
plot((sens.time_myr),(GAST95quant-273),'linewidth',1,'color',c_range)
title('GAST')
xlim([-4e3 0]);
