% plot Imax simulation versus mean field.

% %%% test script 8. Jul
% tmax = 1000;
beta=0.005; gamma=1/40; w=1/100; 
kap=12/10000; delta=2/100;
% betarange = [0.001:0.0001:0.004]; LB = length(betarange);
% betarange = [0.0025]; LB = length(betarange);
% kapparange = [0.0000:0.0002:0.01]; LK = length(kapparange);
% wrange = [0.0000:0.0002:0.01]; LW = length(wrange);
% % rho0=1/100;
% SA = 12; 
% Nt = 1000;
% N = 500;
% mu = 15;
% rho0 = 0.01;
% % 
% % inis = zeros(LB,LW,LK,SA,5);
% % rs = zeros(Nt,LB,LW,LK,SA);
% % is = zeros(Nt,LB,LW,LK,SA);
% imax = zeros(LB,LW,LK,SA);
% rinf = zeros(LB,LW,LK,SA);
% imax_mf = zeros(LB,LW,LK,SA);
% rinf_mf = zeros(LB,LW,LK,SA);
% for wi = 1:LW
%     w = wrange(wi);
%     disp('w')
%     disp(w)
%     for ki = 1:LK
%         kappa = kapparange(ki);
%         disp('kappa')
%         disp(kappa)
%         for bi = 1:LB
%             beta = betarange(bi);
%             disp('beta')
%             disp(beta)
%             
% %             inis_par = ones(SA,5);
%             imax_par = ones(1,SA);
%             rinf_par = ones(1,SA);
%             imax_mf_par = ones(1,SA);
%             rinf_mf_par = ones(1,SA);
%             parfor sa = 1:SA
%                 L=ceil(mu*N/2);
%    
%                 A=ERG(N,L);
%                 x=InitState(rho0,N);
%                 si0=((x==0)*A)*(x==1)'/N;
% 
%                 i0 =sum(x)/N;
%                 rhos=1-i0; rhor=0;
%                 rhoss= (mu/2)-si0;
%                 ini=[rhos,i0,rhor,si0,rhoss];
%                 ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kappa,delta);
%                 [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
%                 result = SIRXi_w_fromAandx(A,x,beta,gamma,w,kappa,delta,tmax,Nt); 
%                 rhoi = result.NIs/result.N;
%                 marhoi = movmean(rhoi,30);
% %                 inis_par(sa) = ini;
% %                 is_par(sa) = result.NIs(1:Nt)/result.N;
% %                 rs_par(sa) = result.NRs(1:Nt)/result.N;
%                 imax_par(sa) = max(marhoi);
%                 rinf_par(sa) = result.NRs(end)/result.N;
%                 imax_mf_par(sa) = max(xs(:,2));
%                 rinf_mf_par(sa) = xs(end,3);
%             end
%             
% %             inis(bi,wi,ki,:,:) = inis_par;
% %             is(1:Nt,bi,wi,ki,:) = is_par;
% %             rs(1:Nt,bi,wi,ki,:) = rs_par;
%             imax(bi,wi,ki,:) = imax_par;
%             rinf(bi,wi,ki,:) = rinf_par;
%             imax_mf(bi,wi,ki,:) = imax_mf_par;
%             rinf_mf(bi,wi,ki,:) = rinf_mf_par;
%         end
%     end
% end


% 
% for beta = betarange
% SA = 5;
% for sa = 1:SA
%     L=ceil(mu*N/2);
%    
%     A=ERG(N,L);
%     x=InitState(rho0,N);
%     si0=((x==0)*A)*(x==1)'/N;
%     
%     i0 =sum(x)/N;
%     rhos=1-i0; rhor=0;
%     rhoss= (mu/2)-si0;
%     ini=[rhos,i0,rhor,si0,rhoss];
%     
%     result{sa} = SIRXi_w_fromAandx(A,x,beta,gamma,w,kap,delta,tmax,Nt); 
%     
% end
% ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);
% [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
% Tt= min(Nt,length(ts));
% 
% diff = 2*xs(:,5)-xs(:,4);
% diffdiff = diff(2:length(diff))-diff(1:(length(diff)-1));
% figure; 
% plot(ts(2:end),diffdiff)
    
% figure;
% for sa = 1:SA
%     rhoi = result{sa}.NIs/result{sa}.N;
%     marhoi = movmean(rhoi,30);
%     pl(sa) = plot(result{sa}.times,marhoi);
%     pl(sa).Color = [0.7,0.7,0.7];
%     pl(sa).LineWidth = 1;
%     hold on;
% end
% pl(SA+1) = plot(ts,xs(:,2));
% pl(SA+1).Color = [0,0,0];
% pl(SA+1).LineWidth = 2;
% pl(SA+1).DisplayName = strcat('rinf is:',num2str(max(xs(:,3))));
% 
% lgd = legend(pl(SA+1));
% lgd.Location = 'northeast';



% load('parametersweep_imax_14Jul.mat')


% SAMPLE RUNS FOR I 
% wi = 2;
% ki = 2;
% bi = 20;
% sa = 5;
% figure;
% for sa=1:SA
%     color = rand(1,3);
%     imaxs = is(:,bi,wi,ki,sa);
%     pl(sa) = plot(imaxs);
%     pl(sa).Color = color;
%     marhoi = movmean(imaxs,30);
%     Imax_ma = max(marhoi) * ones(1,length(imaxs));
%     Imax_mf = imax_mf(bi,wi,ki,sa) * ones(1,length(imaxs));
%     
%     hold on;    
%     pl(sa) = plot(Imax_mf);
%     pl(sa).Color = color;
%     pl(sa)= plot(Imax_ma);
%     pl(sa).Color = color;
%    
% end

% IMAX AGAINST BETA FOR VARIOUS KAPPA AND WS
% take the average over all paths
% 
% for ki=[1,5,10]
%     kappa = kapparange(ki)
%     for wi = [5]
%         w = wrange(wi)
%         imax_mfs_mean = ones(1,length(betarange));
%         imax_sim_mean = ones(1,length(betarange));
%         imax_sim_sdev = ones(1,length(betarange));
%         for bi=1:length(betarange)
%             imax_sim_mean(bi) = nanmean(imax(bi,wi,ki,:));
%             imax_sim_sdev(bi) = nanstd(imax(bi,wi,ki,:));
%             imax_mfs_mean(bi) = nanmean(imax_mf(bi,wi,ki,:));
%         end
%         figure;
%         pl(1) = plot(betarange,imax_mfs_mean);
%         pl(1).Color = 'black';
%         pl(1).LineStyle = '-';
%         pl(1).LineWidth = 2;
%         hold on;
%         pl(2) = errorbar(betarange,imax_sim_mean,imax_sim_sdev);
%         pl(2).Color = 'red';
%         pl(2).LineStyle = 'none';
%         pl(2).LineWidth = 2;
%         pl(2).Marker = 'o';
%         pl(2).MarkerSize = 3;
%         
%         ax = gca;
%         ax.XLabel.String = 'Infection Rate (\beta)';
%         ax.YLabel.String = 'I_{max}';
%         title(strcat('w=',num2str(w),' and kappa=',num2str(kappa)))
%         
%         
%     end
% end





% 
ql = 1;
clear('ql');
bi = 1;
beta = betarange(bi);
wfine=linspace(0.0,0.01,1000);
kapfine=(mu-1)*beta - gamma - wfine;
minval = 0.000;
maxval = 0.08;

LegendFontsizes = 14;
figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
dw=wrange(2)-wrange(1);
dk=kapparange(2)-kapparange(1);
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(kapparange)
        kap = kapparange(ki);
        r = rectangle('Position',[w kap dw dk]');
        
        val = nanmean(imax_mf(bi,wi,ki,:),4);
%         val = atrate_sm(wi,ki);
        if isnan(val)
            r.FaceColor = 'white';
            r.EdgeColor = 'white';
        else
            scaled = (val-minval)/(maxval-minval);
            r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
            r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        end
        r.LineWidth = 0.0001;
    end
end
ql=plot(wfine,kapfine);
ql.LineWidth=3;
ql.Color='black';
ql.DisplayName = 'critical curve (MF)';
ax=gca;
ax.YLim=[0.001,0.01];
ax.XLim=[0.001,0.01];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.25:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.25:1]*(maxval-minval);

lgd=legend(ql);
lgd.Location='northeast';
lgd.Box='off';
lgd.BoxFace.ColorType='truecoloralpha';
lgd.FontSize=LegendFontsizes;

resolution=300;
folder='figures/';
filename=strcat('Simulation_2d_imax_atrate_19Jul');

% direction=strcat(folder,filename,'.png');
% saveas(gcf,direction)




ql = 1;
clear('ql');
bi = 1;
beta = betarange(bi);
wfine=linspace(0.0,0.01,1000);
kapfine=(mu-1)*beta - gamma - wfine;
minval = 0.000;
maxval = 0.08;

LegendFontsizes = 14;
figure; 
plot(wrange,zeros(1,length(wrange)));
hold on;
dw=wrange(2)-wrange(1);
dk=kapparange(2)-kapparange(1);
for wi=1:length(wrange)
    w = wrange(wi);
    for ki=1:length(kapparange)
        kap = kapparange(ki);
        r = rectangle('Position',[w kap dw dk]');
        
        val = nanmean(imax(bi,wi,ki,:),4);
%         val = atrate_sm(wi,ki);
        if isnan(val)
            r.FaceColor = 'white';
            r.EdgeColor = 'white';
        else
            scaled = (val-minval)/(maxval-minval);
            r.FaceColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
            r.EdgeColor = [1,0,0]*scaled +[0,0.5,1]*(1-scaled);
        end
        r.LineWidth = 0.0001;
    end
end
ql=plot(wfine,kapfine);
ql.LineWidth=3;
ql.Color='black';
ql.DisplayName = 'critical curve (MF)';
ax=gca;
ax.YLim=[0.001,0.01];
ax.XLim=[0.001,0.01];
ax.XLabel.String='w';
ax.YLabel.String='\kappa';
colorMap = [linspace(0,1,256)',linspace(0.5,0,256)',linspace(1,0,256)'];
colormap(colorMap);
cb=colorbar;
cb.LineWidth=0.5;
cb.Ticks = [0:0.25:1];
cb.TicksMode ='manual';
cb.TickLabels = minval + [0:0.25:1]*(maxval-minval);

lgd=legend(ql);
lgd.Location='northeast';
lgd.Box='off';
lgd.BoxFace.ColorType='truecoloralpha';
lgd.FontSize=LegendFontsizes;

resolution=300;
folder='figures/';
filename=strcat('Simulation_2d_imax_atrate_19Jul');

% direction=strcat(folder,filename,'.png');
% saveas(gcf,direction)

