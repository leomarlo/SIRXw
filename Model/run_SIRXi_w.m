tmax=2500;
Nt=tmax*2;

N=100;
mu=15;
beta=1/150;
gamma=1/40;  % recovery for the original paper
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

Nrange=[500,1000,1500]; % range for the original paper
wrange=[0:2/10000:1/100]; % range for the original paper
kaprange=[0:2/10000:1/100]; % range for the original paper
brange=[0:0.0002:0.005]; % range for the original paper


Nrange=[500]; % range for the appeal

gamma=1/14;  % recovery for the appeal
wrange=[0:0.01:0.2]; % range for the appeal
kaprange=[0:0.01:0.2]; % range for the appeal
brange=[0:0.004:0.1]; % transmission rate range for the appeal
new_first_wi = 7;

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
    for wi=new_first_wi:length(wrange)
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
            save(strcat('sweepAprl19N',num2str(N),'delta002.mat'),'mfaRinf','mfaImax','simRinf','simImax','wrange','kaprange','brange','N','mu','delta','Nt','tmax')
        end
    end
end
%
% N=100;
% mu=15;
% beta=1/150;
% gamma=1/40;
% w=1/100;
% rho0=1/100;
% kap=12/10000;
% delta=1/100;
% % rhos,rhoi,rhor,rhosi,rhoss
% rhoi=0.01;
% rhos=1-rhoi;
% rhor=0;
% rhosi= rhoi*mu;
% rhoss= (mu/2)-rho0*mu;
% ini=[rhos,rhoi,rhor,rhosi,rhoss];
%
% Nrange=[500,1000,1500];
% wrange=[0:2/10000:1/100];
% kaprange=[0:2/10000:1/100];
% brange=[0:0.0002:0.005];
%
% mfaRinf='to';
% mfaImax='be';
% simRinf='cleared';
% simImax='into nothingness';
%
% Sa=20;
% for Ni=1:length(Nrange)
%     N=Nrange(Ni)
%     clear('mfaRinf')
%     clear('mfaImax')
%     clear('simRinf')
%     clear('simImax')
%     for wi=1:length(wrange)
%         w=wrange(wi)
%         for ki=1:length(kaprange)
%             kap=kaprange(ki)
%             for bi=1:length(brange)
%                 beta=brange(bi);
%                 Imax=zeros(1,Sa);
%                 Rinfs=zeros(1,Sa);
%                 parfor sa=1:Sa
%                     result = SIRXi_w(N,mu,beta,gamma,w,kap,delta,rho0,tmax,Nt);
%     %                 colR{wi,ki,bi,sa}=result.NRs;
%     %                 colX{wi,ki,bi,sa}=result.NXis;
%                     Imax(sa)=max(result.NIs)/result.N;
%                     Rinfs(sa)=max(result.NRs+result.NXis)/result.N;
%                 end
%                 ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);
%                 [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
%                 mfaRinf{wi,ki,bi}=max(1-xs(:,1)-xs(:,2));
%                 mfaImax{wi,ki,bi}=max(xs(:,2));
%                 simRinf{wi,ki,bi}=Rinfs;
%                 simImax{wi,ki,bi}=Imax;
%             end
%             save(strcat('sweepMay24N',num2str(N),'.mat'),'mfaRinf','mfaImax','simRinf','simImax','wrange','kaprange','brange','N','mu','delta','Nt','tmax')
%         end
%     end
% end

load('sweepMay24N200.mat')

wrange=wrange(wrange<0.0044);

mfImax=zeros(length(wrange),length(kaprange),length(brange));
mfRinf=zeros(length(wrange),length(kaprange),length(brange));
meanRinf=zeros(length(wrange),length(kaprange),length(brange));
sdevRinf=zeros(length(wrange),length(kaprange),length(brange));
meanImax=zeros(length(wrange),length(kaprange),length(brange));
sdevImax=zeros(length(wrange),length(kaprange),length(brange));
critical=zeros(length(wrange),length(kaprange));
critical2=zeros(length(wrange),length(kaprange));
critical3=zeros(length(wrange),length(kaprange));

mfImaxatbeta=zeros(length(wrange),length(kaprange));
mfRinfatbeta=zeros(length(wrange),length(kaprange));

nmImaxatbeta=zeros(length(wrange),length(kaprange));
nmRinfatbeta=zeros(length(wrange),length(kaprange));

for wi=1:length(wrange)
    w=wrange(wi)
    for ki=1:length(kaprange)
        kap=kaprange(ki);
        for bi=1:length(brange)
            meanRinf(wi,ki,bi)=mean(simRinf{wi,ki,bi});
            sdevRinf(wi,ki,bi)=std(simRinf{wi,ki,bi});
            meanImax(wi,ki,bi)=mean(simImax{wi,ki,bi});
            sdevImax(wi,ki,bi)=std(simImax{wi,ki,bi});
            mfRinf(wi,ki,bi)=mfaRinf{wi,ki,bi};
            mfImax(wi,ki,bi)=mfaImax{wi,ki,bi};
        end
%         critical(wi,ki)=brange(sdevRinf(wi,ki,:)==max(sdevRinf(wi,ki,:)));
%         critical2(wi,ki)=brange(find(meanRinf(wi,ki,:)>0.1,1,'first'));
%         critical3(wi,ki)=meanRinf(wi,ki,find(brange>=0.002,1,'first'));

        mfImaxatbeta(wi,ki)=mfImax(wi,ki,22);
        mfRinfatbeta(wi,ki)=mfRinf(wi,ki,22);

        nm=simRinf{wi,ki,22};
        nmIm=simImax{wi,ki,22};

        if sum(nm>0.05)>length(nm)/2
            mnm=mean(nm(nm>0.05));
            mnmIm=mean(nmIm(nm>0.05));
        else
            mnm=mean(nm);
            mnmIm=mean(nmIm);
        end
        nmImaxatbeta(wi,ki)=mnmIm;
        nmRinfatbeta(wi,ki)=mnm;

        if mod(ki,10)==0 && mod(wi,10)==0
            figure;
            mfs=zeros(1,length(brange));
            for bi=1:length(brange)
                mfs(bi)=mfaRinf{wi,ki,bi};
            end
            pl=plot(brange,mfs);
            pl.LineWidth=2.5;
            hold on;
            SA=length(simRinf{wi,ki,bi});
            simus=zeros(SA,length(brange));
            for sa=1:SA
                for bi=1:length(brange)
                    temp=simRinf{wi,ki,bi};
                    simus(sa,bi)=temp(sa);
                end
                ql=plot(brange, simus(sa,:));% ql.LineWidth=1;
                ql.LineStyle='None'; ql.Marker='o'; ql.MarkerSize=3; ql.LineWidth=3;
            end
        end

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
image(mfImaxatbeta*200)
ax=gca;
ax.XLabel.String='\kappa';
ax.YLabel.String='w';
ax.XTickLabel=kaplab;
ax.YTickLabel=wlab;

figure;
image(nmImaxatbeta*200)
ax=gca;
ax.XLabel.String='\kappa';
ax.YLabel.String='w';
ax.XTickLabel=kaplab;
ax.YTickLabel=wlab;

figure;
image(mfRinfatbeta*50)
ax=gca;
ax.XLabel.String='\kappa';
ax.YLabel.String='w';
ax.XTickLabel=kaplab;
ax.YTickLabel=wlab;

figure;
image(nmRinfatbeta*50)
ax=gca;
ax.XLabel.String='\kappa';
ax.YLabel.String='w';
ax.XTickLabel=kaplab;
ax.YTickLabel=wlab;

% figure;
% image(critical*5000)
% ax=gca;
% ax.XLabel.String='\kappa';
% ax.YLabel.String='w';
% ax.XTickLabel=kaplab;
% ax.YTickLabel=wlab;

% figure;
% image(critical2*15000)
% ax=gca;
% ax.XLabel.String='\kappa';
% ax.YLabel.String='w';
% ax.XTickLabel=kaplab;
% ax.YTickLabel=wlab;
% % ax.XTicks=
%
% figure;
% image(critical3*1000)
% ax=gca;
% ax.XLabel.String='\kappa';
% ax.YLabel.String='w';
% ax.XTickLabel=kaplab;
% ax.YTickLabel=wlab;

% save('parameter24May.mat','colR','colX','colmaxR','wrange','kaprange','brange','N','mu','delta','Nt','tmax')

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
