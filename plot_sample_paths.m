% clear all 
% 
% load('As_example_1500.mat')
plot_utilities


tmax = 800;
tspan = [0:1:tmax];
Nt = tmax;
SA = 10;

beta = 0.005;
gamma = 1/40;
delta = 1/100;
% rho0 = 0.01;
rhosi0 = 0.15;
kappaws = [[0, 0.0025];[0.0025, 0]; [0.0025,0.0025];[0,0]];

rhoss0 = (mu/2)- rhosi0;
ini=[1-rho0,rho0,0,rhosi0,rhoss0];
% 
% 
% storeI_sm = zeros(tmax+1,length(kappaws), SA);
% storeX_sm = zeros(tmax+1,length(kappaws), SA);
% storeR_sm = zeros(tmax+1,length(kappaws), SA);
% storeSI_sm = zeros(tmax+1,length(kappaws),SA);
% storeI_mf = zeros(tmax+1,length(kappaws));
% storeR_mf = zeros(tmax+1,length(kappaws));
% storeX_mf = zeros(tmax+1,length(kappaws));
% storeSI_mf = zeros(tmax+1,length(kappaws));
% for kwi = 1:length(kappaws)
%     kappa = kappaws(kwi,1);
%     w =  kappaws(kwi,2);
%     disp([kappa,w])
%     ODE = @(y) PA_SIRXi_w(y,beta,gamma,w,kappa,delta);
%     [ts,ys] = ode45(@(t,y) ODE(y),tspan,ini);
%     
%     storeI_mf(:,kwi) = ys(:,2);
%     storeR_mf(:,kwi) = ys(:,3);
%     storeX_mf(:,kwi) = 1 - ys(:,1) - ys(:,2) - ys(:,3);
%     storeSI_mf(:,kwi) = ys(:,4);
%     for sa = 1:SA
%         disp(sa);
%         result=SIRXi_w_fromAandx(A,x',beta,gamma,w,kappa,delta,tmax,Nt);
%         storeI_sm(:,kwi,sa) = result.NIs/result.N;
%         storeR_sm(:,kwi,sa) = result.NRs/result.N;
%         storeX_sm(:,kwi,sa) = result.NXis/result.N;
%         storeSI_sm(:,kwi,sa) = result.NSIs/result.N;
%     end
% end



load('exampleplots18Oct.mat')


% cols = [[0,0,0]; [1,0,0];[0,0,1];[0,1,0]];
for kwi = 1:length(kappaws)
    kappa = kappaws(kwi,1);
    w =  kappaws(kwi,2);
    figure;
    for sa = 1:SA
        pl(sa) = plot(tspan,storeI_sm(:,kwi,sa));
        pl(sa).Color = [0,0,0];
        pl(sa).Color(4) = 0.25;
        pl(sa).LineWidth = 1;
        pl(sa).DisplayName = strcat('Simulation');
        hold on;
    end
    ql(kwi)=plot(tspan,storeI_mf(:,kwi));
    ql(kwi).Color = [0,0,0];
    ql(kwi).LineWidth = 3;
    ql(kwi).DisplayName = strcat('Pair Approximation');
    
    ax = gca;
    ax.XLabel.String = 'time';
    ax.YLabel.String = '\rho_I';
    ax.YLim = [0,0.3];
    ax.XLim = [0,350];    
    ax.YTick = [0:0.1:0.3];
    ax.FontSize = fontsize;
    ax.LineWidth = FrameThickness;
    ax.TickLength=TickLengths;
    ax.FontWeight = FontWeights;

    lgd = legend([ql(kwi),pl(1)]);
    lgd.Location = 'northeast';
    lgd.Box = 'off';
    % lgd.BoxFace.ColorType='truecoloralpha';
    lgd.FontSize=LegendFontsize;

    resolution=300;
    folder='figures/';
    filename=strcat('sampleruns_I_','w_',num2str(w*10000),'_k_',num2str(kappa * 10000),'_', date);
    direction=strcat(folder,filename,'.png');
    saveas(gcf,direction)
end


for kwi = 1:length(kappaws)
    kappa = kappaws(kwi,1);
    w =  kappaws(kwi,2);
    figure;
    for sa = 1:SA
        pl(sa) = plot(tspan,storeR_sm(:,kwi,sa));
        pl(sa).Color = [0,0,0];
        pl(sa).Color(4) = 0.25;
        pl(sa).LineWidth = 1;
        pl(sa).DisplayName = strcat('Simulation');
        hold on;
    end
    ql(kwi)=plot(tspan,storeR_mf(:,kwi));
    ql(kwi).Color = [0,0,0];
    ql(kwi).LineWidth = 3;
    ql(kwi).DisplayName = strcat('Pair Approximation');
    
    ax = gca;
    ax.XLabel.String = 'time';
    ax.YLabel.String = '\rho_R';
    ax.YLim = [0,1];
    ax.XLim = [0,350];
    ax.YTick = [0:0.2:1];
    ax.FontSize = fontsize;
    ax.LineWidth = FrameThickness;
    ax.TickLength=TickLengths;
    ax.FontWeight = FontWeights;

    lgd = legend([ql(kwi),pl(1)]);
    lgd.Location = 'northwest';
    lgd.Box = 'off';
    % lgd.BoxFace.ColorType='truecoloralpha';
    lgd.FontSize=LegendFontsize;

    resolution=300;
    folder='figures/';
    filename=strcat('sampleruns_R_','w_',num2str(w*10000),'_k_',num2str(kappa * 10000),'_', date);
    direction=strcat(folder,filename,'.png');
    saveas(gcf,direction)
end