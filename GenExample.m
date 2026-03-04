for k=1:1000
A=randn(50);
thres=1/2;
A=A-A';
B=expm(A);
P=rand(50)>thres;
C=P.*B;
[U,S,V]=svd(C);
for i=1:10
G(:,:,i)=P.*(U(:,end-i+1)*V(:,end-i+1)');
epsval(i)=S(end-i+1,end-i+1)/norm(G(:,:,i),'fro');
end
[epsort,indsort]=sort(epsval);
if indsort(1)>1
    k
    break
end
end;