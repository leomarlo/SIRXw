function result=SIRXw(N,mu,beta,gamma,w,kap0,kap,rho0,tmax,Nt)
%%%
% takes:
%       infectfract - fraction of initially infected individuals
%       N - number of individuals
%       L - number of edges
%       K - number of runs
%       r - recovery rate
%       l - infection rate
%       w - rewiring rate
%       tmax - maximal  time
%       Nt - number of time data points
%       verbose - enabled if I want to see the timings
%       
% returns:
%       rhoI - number of infected nodes as a function of time
%       drhoI - standard deviation of rhoI
%       rhoSI - number of SI links as a function of time
%       drhoSI - standard deviation of rhoSI
%       rhoSS - number of SS links as a function of time
%       drhoSS - standard deviation of rhoSS
%       times - table with times
%%%

    L=ceil(mu*N/2);
   
    A=ERG(N,L);
    x=InitState(rho0,N);

    NI=sum(x==1);
    NS=N-NI;
    NR=0;
    NX=0;
    NSI=(x*A)*(1-x)';
    NSS=((1-x)*A)*(1-x)'/2;
    NII=(x*A)*x'/2;    
    NSR=0;
    NIR=0;
    NRR=0;    
    NXS=0;
    NXI=0;
    NXR=0;
    NXX=0;
    
    
    rates=[gamma,beta,w,kap0,kap0+kap];
    h=[NI,NSI,NSI,NS,NI];
   
    times=0:tmax/Nt:tmax;
    % initialise time
    t=0;
    % initial time index
    i=1;
    % halting flag is enabled when all nodes are healthy, i.e. susceptible or recovered.
    halt=0;

    % array of N ones.
    N1s=ones(1,N);

    NSs=zeros(1,Nt+1);
    NIs=zeros(1,Nt+1);
    NRs=zeros(1,Nt+1);
    
    NXs=zeros(1,Nt+1);
    
    NSSs=zeros(1,Nt+1);
    NSIs=zeros(1,Nt+1);
    NSRs=zeros(1,Nt+1);
    NIIs=zeros(1,Nt+1);
    NIRs=zeros(1,Nt+1);
    NRRs=zeros(1,Nt+1);
    NXSs=zeros(1,Nt+1);
    NXIs=zeros(1,Nt+1);
    NXRs=zeros(1,Nt+1);
    NXXs=zeros(1,Nt+1);

    while t<tmax 
        a=h.*rates;
        as=sum(a);
        r1=rand(); r2=rand();
        acum=a(1);
        cond1=r1*as>=acum; acum=acum+a(2);
        cond2=r1*as>=acum; acum=acum+a(3);
        cond3=r1*as>=acum; acum=acum+a(4);
        cond4=r1*as>=acum;
        reactiontype = 1+cond1+cond2+cond3+cond4;
        tau=log(1/r2)/as;
        
        switch reactiontype
            case 1
                % recovery
                nI=find(x==1);
        %         disp(?nI)
                if isempty(nI)
                    tau=0;
                else
                    tbr=randsample([nI,nI],1); % to be recovered
                    x(tbr)=2;
                end

            case 2
                % infection
                SI=(A & ((x==0)')*N1s & (N1s')*(x==1));
                % indices of SI-links
                [rSI,cSI]=find(SI);

                x(randsample([rSI;rSI],1))=1;

            case 3
                % rewiring

                SI=(A & ((x==0)')*N1s & (N1s')*(x==1));
                % indices of SI-links
                [rSI,cSI]=find(SI);
                rn=randi(length(rSI));

                examples=find(A(rSI(rn),:)==0 & x~=1 & (1:N)~=rSI(rn));
                if isempty(examples)
                    % t shouldnt go forward
                    tau=0;
                else
                    re=randsample([examples,examples],1);
                    A(rSI(rn),cSI(rn))=0;
                    A(cSI(rn),rSI(rn))=0;

                    A(rSI(rn),re)=1;
                    A(re,rSI(rn))=1;

                end

            case 4
                % Susceptible is quarantined
                nS=find(x==0);
        %         disp(?nI)
                if isempty(nS)
                    tau=0;
                else
                    tbr=randsample([nS,nS],1); % to be recovered
                    x(tbr)=3; % TODO CHANGE TO ANOTHER VALUE
                end
                    

            case 5
                % Infected is quarantied
                nI=find(x==1);
        %         disp(?nI)
                if isempty(nI)
                    tau=0;
                else
                    tbr=randsample([nI,nI],1); % to be recovered
                    x(tbr)=3; % TODO CHANGE TO ANOTHER VALUE
                end
            otherwise
                warning('cannot be')
                
        end

        % recalculate A,SS,SI,II,SR,IR,RR, I, S and R
        NSS=((x==0)*A)*(x==0)'/2;
        NSI=((x==0)*A)*(x==1)';
        NSR=((x==0)*A)*(x==2)';
        NII=((x==1)*A)*(x==1)'/2;
        NIR=((x==1)*A)*(x==2)';
        NRR=((x==2)*A)*(x==2)'/2;

        NXS=((x==3)*A)*(x==0)';
        NXI=((x==3)*A)*(x==1)';
        NXR=((x==3)*A)*(x==2)';
        NXX=((x==3)*A)*(x==3)'/2;

        NS=sum(x==0);
        NI=sum(x==1);
        NR=sum(x==2);
        NX=sum(x==3);

        % update h
        h=[NI,NSI,NSI,NS,NI];
        % update time
        t=t+tau;

        ind=find(times<=t,1,'last');

        % if ind>i
        if ind>=i+1
            for k=[i+1:ind]
                % update the observables 
                % recalculate A,SS,SI,II,SR,IR,RR, I, S and R

                NSs(k)=NS;
                NIs(k)=NI;
                NRs(k)=NR;
                
                NXs(k)=NX;
                NSSs(k)=NSS;
                NSIs(k)=NSI;
                NSRs(k)=NSR;
                NIIs(k)=NII;
                NIRs(k)=NIR;
                NRRs(k)=NRR;
                NXSs(k)=NXS;
                NXIs(k)=NXI;
                NXRs(k)=NXR;
                NXXs(k)=NXX;
            end    
        end

        %update the new index.
        i=ind;  
        if NI==0
            halt=1;
            for k=[i+1:Nt+1]
                NSs(k)=NS;
                NIs(k)=NI;
                NRs(k)=NR;
                NXs(k)=NX;
                NSSs(k)=NSS;
                NSIs(k)=NSI;
                NSRs(k)=NSR;
                NIIs(k)=NII;
                NIRs(k)=NIR;
                NRRs(k)=NRR;
                NXSs(k)=NXS;
                NXIs(k)=NXI;
                NXRs(k)=NXR;
                NXXs(k)=NXX;
            end
            break
        end

    end


    result.NSs=NSs;
    result.NIs=NIs;
    result.NRs=NRs;
    result.NXs=NXs;
    result.NSSs=NSSs;
    result.NSIs=NSIs;
    result.NSRs=NSRs;
    result.NIIs=NIIs;
    result.NIRs=NIRs;
    result.NRRs=NRRs;
    result.NXSs=NXSs;
    result.NXIs=NXIs;
    result.NXRs=NXRs;
    result.NXXs=NXXs;
    result.N=N;
    result.L=L;
    result.beta=beta;
    result.gamma=gamma;
    result.w=w;
    result.mu=mu;
    result.kap0=kap0;
    result.kap=kap;
    result.tmax=tmax;
    result.Nt=Nt;
    result.halt=halt;
    result.times=times;
    result.rates=rates;
    
end