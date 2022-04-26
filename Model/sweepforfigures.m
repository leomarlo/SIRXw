%% plot 20 runs above the treshold
%% small kappa, large kappa, positive delta, zero delta, small w, large w, N = 500, mu = 15, mu = 5,

tmax=500;
Nt=tmax*2;
N=1000;
mu=15;
L = ceil(mu*N/2);
beta=0.004;
gamma=1/40;
w=1/100;
rho0=1/100;
kap=12/10000;
delta=2/100;
Ns = [1000,2000];
kappas = [0.001,0.01];
ws = [0.001,0.01];
deltas = [4/100];
i0s = [2/500,10/500];
SA = 3;
same = 8;
store_sm = zeros(length(Ns),length(kappas),length(ws),length(i0s),SA,same,15,Nt+1);
store_mf = zeros(length(Ns),length(kappas),length(ws),length(i0s),SA,6,Nt+1);

Ni=1;
for N=Ns
    ki=1;
    for kappa=kappas
        wi=1;
        for w=ws
            di=1;
            for delta=deltas
                i0i=1;
                for i0=i0s
                    for sa=1:SA
                        A=ERG(N,L);
                        x=InitState(i0,N);

                        rhos=1-i0; rhor=0;
                        si0=((x==0)*A)*(x==1)'/N;
                        rhoss= (mu/2)-si0;
                        ini=[rhos,i0,rhor,si0,rhoss];
                        ODE = @(x) PA_SIRXi_w(x,beta,gamma,w,kap,delta);
                        [ts,xs] = ode45(@(t,x) ODE(x),[0 tmax],ini);
                        Tt= min(Nt,length(ts));
                        store_mf(Ni,ki,wi,i0i,sa,1,1:Tt)=ts(1:Tt);
                        store_mf(Ni,ki,wi,i0i,sa,2,1:Tt)=xs(1:Tt,1); %s   % ssi = 2ss si /s = 2 * [4] *[5] /1
                        store_mf(Ni,ki,wi,i0i,sa,3,1:Tt)=xs(1:Tt,2); %i
                        store_mf(Ni,ki,wi,i0i,sa,4,1:Tt)=xs(1:Tt,3); %r
                        store_mf(Ni,ki,wi,i0i,sa,5,1:Tt)=xs(1:Tt,4); %si
                        store_mf(Ni,ki,wi,i0i,sa,6,1:Tt)=xs(1:Tt,5); %ss
                        disp(strcat('kap=',num2str(kappa),' w=',num2str(w),' delta=',num2str(delta),' i0=',num2str(i0),' si0=',num2str(si0),' sa=',num2str(sa)));
                        for sm=1:same

                            result = SIRXi_w_fromAandx(A,x,beta,gamma,w,kap,delta,tmax,Nt); 
    %                         store_sm(ki,wi,ws,di,i0i,sa,sm,1,:)=
                            store_sm(Ni,ki,wi,i0i,sa,sm,1,:)=result.times;
                            store_sm(Ni,ki,wi,i0i,sa,sm,2,:)=result.NSs; %1
                            store_sm(Ni,ki,wi,i0i,sa,sm,3,:)=result.NIs; %2
                            store_sm(Ni,ki,wi,i0i,sa,sm,4,:)=result.NRs; %3
                            store_sm(Ni,ki,wi,i0i,sa,sm,7,:)=result.NXis; %6
                            store_sm(Ni,ki,wi,i0i,sa,sm,6,:)=result.NSSs; %5
                            store_sm(Ni,ki,wi,i0i,sa,sm,5,:)=result.NSIs; %4
                            store_sm(Ni,ki,wi,i0i,sa,sm,8,:)=result.NSRs; %7
                            store_sm(Ni,ki,wi,i0i,sa,sm,9,:)=result.NIIs; %8
                            store_sm(Ni,ki,wi,i0i,sa,sm,10,:)=result.NIRs; %9
                            store_sm(Ni,ki,wi,i0i,sa,sm,11,:)=result.NRRs; %10
                            store_sm(Ni,ki,wi,i0i,sa,sm,12,:)=result.NXiSs; %11
                            store_sm(Ni,ki,wi,i0i,sa,sm,13,:)=result.NXiIs; %12
                            store_sm(Ni,ki,wi,i0i,sa,sm,14,:)=result.NXiRs; %13
                            store_sm(Ni,ki,wi,i0i,sa,sm,15,:)=result.NXiXis; %14
    %                         checkitout
                            u = 1;
                        end
                        u = 1;
                    end

                    i0i=i0i+1;
                end
                di=di+1;
            end
            wi=wi+1;
        end
        ki=ki+1;
        save('samplepaths23Jun.mat','store_sm','store_mf','beta','gamma','w','kap','delta','N','tmax','Nt','mu','L','Ns','kappas','ws','deltas','i0s') 
    end
    save('samplepaths23Jun.mat','store_sm','store_mf','beta','gamma','w','kap','delta','N','tmax','Nt','mu','L','Ns','kappas','ws','deltas','i0s') 
    Ni=Ni+1;
end

figure;
Ni=2; ki=1; wi=1; i0i=1; sa=3; ind = 3;
for sm=1:same
    times= reshape(store_sm(Ni,ki,wi,i0i,sa,sm,1,:),Nt+1,[]);
    values=reshape(store_sm(Ni,ki,wi,i0i,sa,sm,ind,:),Nt+1,[]);
    plot(times,values/Ns(Ni),'LineWidth',1,'Color','blue') 
    hold on;
end
times_mf = reshape(store_mf(Ni,ki,wi,i0i,sa,1,:),Nt+1,[]);
values_mf = reshape(store_mf(Ni,ki,wi,i0i,sa,ind,:),Nt+1,[]);
LT = find(times_mf>0,1,'last');
plot(times_mf(1:LT),values_mf(1:LT),'LineWidth',3,'Color','black')

