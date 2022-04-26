function result=SIRXi_w(N,mu,beta,gamma,w,kap,delta,rho0,tmax,Nt)



    L=ceil(mu*N/2);
   
    A=ERG(N,L);
    x=InitState(rho0,N);

    NI=sum(x==1); % 1
    NS=N-NI; % 0
    NR=0; % 2
    NXi=0; % 3
    NSS=((x==0)*A)*(x==0)'/2;
    NSI=((x==0)*A)*(x==1)';
    NII=((x==1)*A)*(x==1)'/2;  
    NSR=0; NIR=0; NRR=0; NXiS=0; NXiI=0; NXiR=0; NXiXi=0; 
    
    
    rates=[gamma,beta,w,kap,delta];
    h=[NI,NSI,NSI,NI,NXi];
   
    times=0:tmax/Nt:tmax;
    t=0; % initialise time
    i=1; % initial time index
    halt=0; % halting flag is enabled when all nodes are healthy, i.e. susceptible or recovered.
    
    
    N1s=ones(1,N);

    NSs=zeros(1,Nt+1);
    NIs=zeros(1,Nt+1);
    NRs=zeros(1,Nt+1);
    NXss=zeros(1,Nt+1);
    NXis=zeros(1,Nt+1);
                
    NSSs=zeros(1,Nt+1);
    NSIs=zeros(1,Nt+1);
    NSRs=zeros(1,Nt+1);
    NIIs=zeros(1,Nt+1);
    NIRs=zeros(1,Nt+1);
    NRRs=zeros(1,Nt+1);

    NXiSs=zeros(1,Nt+1);
    NXiIs=zeros(1,Nt+1);
    NXiRs=zeros(1,Nt+1);
    NXiXis=zeros(1,Nt+1);
    
%     cases=[];
    
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
%         cases(end+1)=reactiontype;
        
        switch reactiontype
            case 1
                % recovery
                nI=find(x==1);
                tbr=randsample([nI,nI],1); % to be recovered
                x(tbr)=2;

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
                % Infected is quarantied
                nI=find(x==1);
                tbr=randsample([nI,nI],1); % to be recovered
                x(tbr)=3; % TODO CHANGE TO ANOTHER VALUE
             case 5
                % recovered infected is resocialized
                nXi=find(x==3);
                x(randsample([nXi,nXi],1))=0;
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

        NXiS=((x==3)*A)*(x==0)';
        NXiI=((x==3)*A)*(x==1)';
        NXiR=((x==3)*A)*(x==2)';
        NXiXi=((x==3)*A)*(x==3)'/2;

        NS=sum(x==0);
        NI=sum(x==1);
        NR=sum(x==2);
        NXi=sum(x==3);

        % update h
        h=[NI,NSI,NSI,NI,NXi];
        % update time
        t=t+tau;

        ind=find(times<=t,1,'last');
    
        if ind>=i+1
            for k=[i+1:ind]
                % update the observables 
                % recalculate A,SS,SI,II,SR,IR,RR, I, S and R

                NSs(k)=NS;
                NIs(k)=NI;
                NRs(k)=NR;
                NXis(k)=NXi;
                NSSs(k)=NSS;
                NSIs(k)=NSI;
                NSRs(k)=NSR;
                NIIs(k)=NII;
                NIRs(k)=NIR;
                NRRs(k)=NRR;
                
                NXiSs(k)=NXiS;
                NXiIs(k)=NXiI;
                NXiRs(k)=NXiR;
                NXiXis(k)=NXiXi;
                
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
                NXis(k)=NXi;
                NSSs(k)=NSS;
                NSIs(k)=NSI;
                NSRs(k)=NSR;
                NIIs(k)=NII;
                NIRs(k)=NIR;
                NRRs(k)=NRR;
                
                NXiSs(k)=NXiS;
                NXiIs(k)=NXiI;
                NXiRs(k)=NXiR;
                NXiXis(k)=NXiXi;
            end
            break
        end

    end
    
    
    result.NSs=NSs;
    result.NIs=NIs;
    result.NRs=NRs;
    result.NXis=NXis;
    result.NSSs=NSSs;
    result.NSIs=NSIs;
    result.NSRs=NSRs;
    result.NIIs=NIIs;
    result.NIRs=NIRs;
    result.NRRs=NRRs;

    result.NXiSs=NXiSs; % (k)=NXiS;
    result.NXiIs=NXiIs;
    result.NXiRs=NXiRs;
    result.NXiXis=NXiXis;
    result.N=N;
    result.L=L;
    result.beta=beta;
    result.gamma=gamma;
    result.w=w;
    result.delta=delta;
    result.mu=mu;
    result.kap=kap;
    result.tmax=tmax;
    result.Nt=Nt;
    result.halt=halt;
    result.times=times;
    result.rates=rates;
%     result.cases=cases;
    
end
