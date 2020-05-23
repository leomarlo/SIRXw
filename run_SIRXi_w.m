tmax=500;
Nt=1000;

N=700;
mu=14; 
beta=1/150;
gamma=1/40;
w=1/100;
rho0=1/100;
kap=12/10000;
delta=1/100;
% rhos,rhoi,rhor,rhosi,rhoss
rhoi=0.014;
rhos=1-rhoi;
rhor=0;
rhosi= rhoi*mu;
rhoss= (mu/2)-rho0*mu;
ini=[rhos,rhoi,rhor,rhosi,rhoss];

[beta,(w+gamma+kap)/(mu+1)]

(gamma+2/100)/(mu+1)
wrange=[0:2/1000:1/100];
kaprange=[0:2/1000:1/100];
brange=[0:0.0005:0.006];

for wi=1:length(wrange)
    w=wrange(wi)
    for ki=1:length(kaprange)
        kap=kaprange(ki)
        for bi=1:length(brange)
            beta=brange(bi)
            for sa=1:20
                result = SIRXi_w(N,mu,beta,gamma,w,kap,delta,rho0,tmax,Nt);
                colR{wi,ki,bi,sa}=result.NRs;
                colX{wi,ki,bi,sa}=result.NXis;
                colmaxR{wi,ki,bi,sa}=max(result.NRs);
                colmaxI{wi,ki,bi,sa}=max(result.NIs);
                
            end
            ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);
            [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
            
        end
    end
end


meanmax=zeros(length(wrange),length(kaprange),length(brange));
sdevmax=zeros(length(wrange),length(kaprange),length(brange));
critical=zeros(length(wrange),length(kaprange));
critical2=zeros(length(wrange),length(kaprange));

for wi=1:length(wrange)
    w=wrange(wi);
    for ki=1:length(kaprange)
        kap=kaprange(ki);
        for bi=1:length(brange)
            vec=[colmaxR{wi,ki,bi,1},colmaxR{wi,ki,bi,2},colmaxR{wi,ki,bi,3},colmaxR{wi,ki,bi,4}]/N;
            meanmax(wi,ki,bi)=mean(vec);
            sdevmax(wi,ki,bi)=std(vec);
        end
        critical(wi,ki)=brange(sdevmax(wi,ki,:)==max(sdevmax(wi,ki,:)));
        critical2(wi,ki)=brange(find(meanmax(wi,ki,:)>0.1,1,'first'));
    end
end

kaplab={};
wlab={};
for ki=1:length(kaprange)
    kap=kaprange(ki);
    kaplab{ki}=num2str(kap);
end
for wi=1:length(wrange)
    w=wrange(wi);
    wlab{wi}=num2str(w);
end

figure;
image(critical2*15000)
ax=gca;
ax.XLabel.String='\kappa';
ax.YLabel.String='w';
ax.XTickLabel=kaplab;
ax.YTickLabel=wlab;

% save('parametersweep23May.mat','colR','colX','colmaxR','wrange','kaprange','brange','N','mu','delta','Nt','tmax')

% figure
% for i=1:5
%     disp('he')
% %     disp(i)
%     result{i} = SIRXi_w(N,mu,beta,gamma,w,kap,delta,rho0,tmax,Nt);
% end  
% ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);

% %% ODE 45 integration routine.
% [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
% 
% resolution=300;
% folder='figures/';
% %% Compare to MODEL 1
% filename=strcat('22May_MF_vs_simulation');
% figure;
% pl(1)=plot(ts,xs(:,2));
% pl(1).LineWidth=3;
% pl(1).DisplayName='Mean Field';
% hold on;
% for i=1:5
%     pl(i+1)=plot(result{i}.times, result{i}.NIs/result{i}.N);
%     pl(i+1).DisplayName=strcat('Simulation ',num2str(i));
% end
% lgd=legend(pl);

% direction=strcat(folder,filename,'.png');
% saveas(gcf,direction)
% print(direction,'-dpng',resolution);