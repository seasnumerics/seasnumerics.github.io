function F=Force(X)
global nv K l1 l2;

F = zeros(nv,3);
lmax = nv + 5120-2;
for l=1:lmax
    Flink=K*(X(l2(l),:)-X(l1(l),:));
    F(l1(l),:)=F(l1(l),:)+Flink;
    F(l2(l),:)=F(l2(l),:)-Flink;
end