function [y,a]=plane_project2(x,E)
    [N,M]=size(x);
    p=size(E,2);
    a=zeros(p,M);
    ct=E(:,1);
    Ep=E(:,2:p)-ct*ones(1,p-1);
    a(2:p,:)=Ep\(x-ct*ones(1,M));
    a(1,:)=ones(1,M)-sum(a(2:p,:),1);
    y=E*a;
end