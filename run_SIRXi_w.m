tmax=2500;
Nt=tmax*2;

N=100;
mu=15; 
beta=1/150;
gamma=1/40;
w=1/100;
rho0=1/100;
kap=12/10000;
delta=2/100;
% rhos,rhoi,rhor,rhosi,rhoss
rhoi=0.01;
rhos=1-rhoi;
rhor=0;
rhosi= rhoi*mu;
rhoss= (mu/2)-rho0*mu;
ini=[rhos,rhoi,rhor,rhosi,rhoss];

Nrange=[500,1000,1500];
wrange=[0:2/10000:1/100];
kaprange=[0:2/10000:1/100];
brange=[0:0.0002:0.005];

mfaRinf='to';
mfaImax='be';
simRinf='cleared';
simImax='into nothingness';

Sa=20;
for Ni=1:length(Nrange)
    N=Nrange(Ni)
    clear('mfaRinf')
    clear('mfaImax')
    clear('simRinf')
    clear('simImax')
    for wi=1:length(wrange)
        w=wrange(wi)
        for ki=1:length(kaprange)
            kap=kaprange(ki)
            for bi=1:length(brange)
                beta=brange(bi);
                Imax=zeros(1,Sa);
                Rinfs=zeros(1,Sa);
                parfor sa=1:Sa
                    result = SIRXi_w(N,mu,beta,gamma,w,kap,delta,rho0,tmax,Nt);
    %                 colR{wi,ki,bi,sa}=result.NRs;
    %                 colX{wi,ki,bi,sa}=result.NXis;
                    Imax(sa)=max(result.NIs)/result.N;
                    Rinfs(sa)=max(result.NRs+result.NXis)/result.N;
                end
                ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);
                [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
                mfaRinf{wi,ki,bi}=max(1-xs(:,1)-xs(:,2));
                mfaImax{wi,ki,bi}=max(xs(:,2));
                simRinf{wi,ki,bi}=Rinfs;
                simImax{wi,ki,bi}=Imax;
            end
            save(strcat('sweepJun06N',num2str(N),'delta002.mat'),'mfaRinf','mfaImax','simRinf','simImax','wrange','kaprange','brange','N','mu','delta','Nt','tmax')
        end
    end
end
% 
% meanmax=zeros(length(wrange),length(kaprange),length(brange));
% sdevmax=zeros(length(wrange),length(kaprange),length(brange));
% critical=zeros(length(wrange),length(kaprange));
% critical2=zeros(length(wrange),length(kaprange));
% 
% for wi=1:length(wrange)
%     w=wrange(wi);
%     for ki=1:length(kaprange)
%         kap=kaprange(ki);
%         for bi=1:length(brange)
%             vec=[colmaxR{wi,ki,bi,1},colmaxR{wi,ki,bi,2},colmaxR{wi,ki,bi,3},colmaxR{wi,ki,bi,4}]/N;
%             meanmax(wi,ki,bi)=mean(vec);
%             sdevmax(wi,ki,bi)=std(vec);
%         end
%         critical(wi,ki)=brange(sdevmax(wi,ki,:)==max(sdevmax(wi,ki,:)));
%         critical2(wi,ki)=brange(find(meanmax(wi,ki,:)>0.1,1,'first'));
%     end
% end
% 
% kaplab={};
% wlab={};
% for ki=1:length(kaprange)
%     kap=kaprange(ki);
%     kaplab{ki}=num2str(kap);
% end
% for wi=1:length(wrange)
%     w=wrange(wi);
%     wlab{wi}=num2str(w);
% end
% 
% figure;
% image(critical2*15000)
% ax=gca;
% ax.XLabel.String='\kappa';
% ax.YLabel.String='w';
% ax.XTickLabel=kaplab;
% ax.YTickLabel=wlab;

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