function x=InitState(ra,N)
K=floor(N*ra); % number of ones
L=N-K; %number of zeros
xold=ones(1,K);
for k=(1:L)
    K=length(xold);
    rd=randi([0,K]);
    xnew=[xold(1:rd),0,xold(rd+1:K)];
    xold=xnew;
end
x=xnew;
end
