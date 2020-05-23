function A=ERG(N,L)
kd=eye(N);
A=sparse(zeros(N));
[rt,c]=find(kd==0);
r=rt(rt<c);
c=c(rt<c);
K=length(r);
edges=[r,c];
choise=randperm(K,L);
chosenedges=edges(choise,:);
linearindices=sub2ind(size(A), chosenedges(:,1),chosenedges(:,2));
A(linearindices)=1;
A=sparse(A+A');

end