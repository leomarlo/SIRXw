plot_utilities

% 
% betarange = 0.00001:0.00004:0.0025;
% LB = length(betarange);
% % beta=0.004; 
% gamma=1/40; 
% w=1/100; 
% kap=12/10000; 
% delta=2/100;
% mu = 15;
% rho0 = 0.01;
% rhosi0 = rho0 * mu;
% rhoss0 = (mu/2) - rhosi0;
% ini=[1-rho0,rho0,0,rhosi0,rhoss0];
% rates = w + kap + gamma;
% 
% load('As.mat')
% sizaAs = size(As);
% 
% estimated = zeros(1,length(betarange));
% mean_field = zeros(1,length(betarange));
% simulations = zeros(length(betarange),sizaAs(3));
% simulation = zeros(length(betarange),sizaAs(3));
% simulation_mean = zeros(1,length(betarange));
% simulation_sdev = zeros(1,length(betarange));
% 
% 
% for bi = 1:LB
%     beta = betarange(bi);
%     sizaAs = size(As);
%     sims = zeros(1,sizaAs(3));
%     parfor sa = 1:sizaAs(3)
%         A = As(:,:,sa);
%         x = xs(:,sa);
%         result = SIRXi_w_fromAandx(A,x',beta,gamma,w,kap,delta,tmax,tmax);
%         rhorinf = result.NRs(end)/result.N;
%         sims(sa) = rhorinf;
%     end
%     simulations(bi,:) = sims;
%     simulation_mean(bi) = nanmean(sims);
%     simulation_sdev(bi) = std(sims);
%     disp(bi);
% end
% 
% for bi = 1:LB
%     beta =  betarange(bi);
%     fun = @(rinf) (rinf - implicitrinf(beta,rates,rho0,rhosi0,rhoss0,rinf));
%     rinf0 = 0.05;
%     rinf_optimum = fzero(fun,rinf0);
%     estimated(bi) = rinf_optimum;
%     ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);
%     [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
%     mean_field(bi) = xs(end,3);
% end

load('plot_implicit_function_02Sept.mat')
sisi = size(simulations);
pl = 1;


clear('pl');
figure; 
pl(3) = errorbar(betarange,simulation_mean, simulation_sdev/sqrt(sisi(2)));
% pl(3) = plot(betarange,simulations);
pl(3).Color = 'blue';
pl(3).LineWidth = 1.75;
pl(3).LineStyle = 'none';
pl(3).Marker = 'o';
pl(3).MarkerSize = 5;
pl(3).DisplayName = 'Simulation';
hold on;
pl(1) = plot(betarange,mean_field);
pl(1).LineWidth = 2.75;
pl(1).Color = 'black';
pl(1).DisplayName = 'Pair approximation';
% hold on;
pl(2) = plot(betarange,estimated);
pl(2).Color = 'red';
pl(2).LineWidth = 2.75;
pl(2).LineStyle = '--';
pl(2).DisplayName = 'Implicit equation from initial cond.';
% 
% pl(2).DisplayName = 'Implicit equation from \rho_S(0) and \rho_{SI}(0)';

ax = gca;
ax.YLabel.String = 'r_\infty';
ax.XLabel.String = '\beta';
ax.YLim = [0,0.12];
ax.XLim = [0,0.00252];
ax.FontSize = fontsize;
ax.LineWidth = FrameThickness;
ax.TickLength=TickLengths;
ax.FontWeight = FontWeights;

% daspect([1 1 1])

lgd = legend(pl);
lgd.Location = 'northwest';
% lgd.BoxFace.ColorType='truecoloralpha';
lgd.FontSize=LegendFontsize;

% grid on
% daspect([1 1 1])


resolution=600;
folder='figures/';
filename=strcat('Rinf_comparison_', date);
direction=strcat(folder,filename,'.png');
saveas(gcf,direction)
% 
% sisi = size(simulations);
% figure;
% db = betarange(2)-betarange(1);
% for jj = 1:sisi(2)
%     plt = plot(betarange + db*0.001*jj,simulations(:,jj));
%     plt.Marker = '.';
%     plt.MarkerSize = 4.5;
%     plt.LineStyle = 'none';
%     hold on;
% end
