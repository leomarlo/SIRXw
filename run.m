tmax=10000;
Nt=1000;

N=1000;
mu=6;
beta=0.001;
gamma=0.002;
w=0.001;
rho0=0.02;

%% COMPARE PA WITH SIMULATION

modl.f1=@(x) PA(x,1,beta,gamma,w);
modl.f2=@(x) PA(x,2,beta,gamma,w);
modl.f20=@(x) PA(x,2,beta,gamma,0);

ini=[rho0,0,rho0*mu,(mu/2)-rho0*mu];

%% ODE 45 integration routine.
[ts1,xs1] = ode45(@(t,x) modl.f1(x),[0 tmax],ini);
[ts2,xs2] = ode45(@(t,x) modl.f2(x),[0 tmax],ini);

result0=SIRw(N,mu,beta,gamma,0,rho0,2,tmax,Nt);
resultw=SIRw(N,mu,beta,gamma,w,rho0,2,tmax,Nt);


%% Simulation
%% say a few simulations
UU=5;
for i=1:UU
    result1(i)=SIRw(N,mu,beta,gamma,w,rho0,1,tmax,Nt);
    result2(i)=SIRw(N,mu,beta,gamma,w,rho0,2,tmax,Nt);
end

% load('simulation_results.mat')

folder='figures/';
%% Compare to MODEL 1
Model=1;
filename=strcat('comparison_08May_',num2str(Model));
LineWidths=1;
% clear('fig');clear('lgd');clear('pl');
fig=figure;
for i=1:UU
    pl(i)=plot(result1(i).times,result1(i).NIs/result1(i).N);
    hold on
    pl(i).LineWidth=LineWidths;
    p1(i).DisplayName='Simulation for ';
end
pl(UU+1)=plot(ts1,xs1(:,1));
% pl(UU+1).LineType=2;
pl(UU+1).LineWidth=LineWidths*3;
pl(UU+1).Color='black';
p1(UU+1).DisplayName='Mean Field For Model 1';
ax=gca;
ax.XLabel.String='Time';
ax.YLabel.String='Prevalence';
lgd=legend(pl);
lgd.Location='northeast';

direction=strcat(folder,filename,'.png');
print(direction,'-dpng',resolution);

%% Compare to MODEL 2
Model=2;
filename=strcat('comparison_08May_',num2str(Model));
LineWidths=1;
% clear('fig');clear('lgd');clear('pl');
fig=figure;
for i=1:UU
    pl(i)=plot(result2(i).times,result2(i).NIs/result2(i).N);
    hold on
    pl(i).LineWidth=LineWidths;
    p1(i).DisplayName='Simulation for ';
end
pl(UU+1)=plot(ts2,xs2(:,1));
% pl(UU+1).LineType=2;
pl(UU+1).LineWidth=LineWidths*3;
pl(UU+1).Color='black';
p1(UU+1).DisplayName='Mean Field For Model 1';
ax=gca;
ax.XLabel.String='Time';
ax.YLabel.String='Prevalence';
lgd=legend(pl);
lgd.Location='northeast';

direction=strcat(folder,filename,'.png');
print(direction,'-dpng',resolution);

% load('Model_',num2str(Model),'_08May.mat'); % loads sweep.







% lgd.Orientation='vertical';
% lgd.Box='off';
% lgd.FontSize=LegendFontsizes*(0.5);
% lgd.Color=blackcolor;
% lgd.LineWidth=LineWidths;

%% RUN A PARAMETER SWEEP

% wrange=0:0.0001:0.01;
% K=10;
% for Model=1:2
%     maxI=zeros(length(wrange),K);
%     limR=zeros(length(wrange),K);
%     for wind=1:length(wrange)
%         w=wrange(wind);
%         disp(w)
%         for sa=1:K
%             disp(sa)
%             result=SIRw(N,mu,beta,gamma,w,rho0,Model,tmax,Nt);
%             maxI(wind,sa)=max(result.NIs/N);
%             limR(wind,sa)=max(result.NRs/N);
% 
%         end
%     end
% 
%     sweep.maxI=maxI;
%     sweep.Rinf=limR;
% 
%     folder='figures/';
%     filename=strcat('Rinf_w_sweep_08May_Model_',num2str(Model));
%     resol=300;
%     resolution=strcat('-r',num2str(resol));
% 
%     xs=sweep.Rinf;
%     fig=figure;
%     er=errorbar(wrange,mean(xs,2),std(xs'));
%     er.Marker='o';
%     er.MarkerSize=6;
%     ax=gca;
%     ax.YLim=[0,1];
%     ax.XLim=[0,max(wrange)];
%     ax.XLabel.String='w';
%     ax.YLabel.String='r_{\infty}';
% 
%     direction=strcat(folder,filename,'.png');
%     print(direction,'-dpng',resolution);
% %     savefig(fig,strcat(direction(1:(end-4)),'.fig'));
% 
%     folder='figures/';
%     filename=strcat('maxI_w_sweep_07May_Model_',num2str(Model));
%     resol=300;
%     resolution=strcat('-r',num2str(resol));
% 
%     xs=sweep.maxI;
%     fig=figure;
%     er=errorbar(wrange,mean(xs,2),std(xs'));
%     er.Marker='o';
%     er.MarkerSize=6;
%     ax=gca;
%     ax.YLim=[0,0.5];
%     ax.XLim=[0,max(wrange)];
%     ax.XLabel.String='w';
%     ax.YLabel.String='max(\rho_I)';
% 
%     direction=strcat(folder,filename,'.png');
%     print(direction,'-dpng',resolution);
% %     savefig(fig,strcat(direction(1:(end-4)),'.fig'));
% 
%     save(strcat('Model_',num2str(Model),'_09May.mat'),'sweep')
% end

