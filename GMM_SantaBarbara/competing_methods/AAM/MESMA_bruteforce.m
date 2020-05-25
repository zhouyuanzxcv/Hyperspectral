function [idx, A, rec, minerr]=MESMA_bruteforce (x,L)
% MESMA_BRUTEFORCE Efficient brute force approach for MESMA problems with
% a low number of endmember libraries.
%
% input: x contains the mixed spectra (dimension x number of spectra)
%        L contains the libraries as a cell array containing p matrices
%          each matrix is of size (dimension x whatever)
% output: idx contains the indices into the libraries of the best endmembers
%         A is the abundance matrix wrt. these endmembers
%         rec contains the spectral reconstructions
%         minerr is the reconstruction error (RMS error)
%


[d,num]=size(x);
p=numel(L);
for i=1:p
    N(i)=size(L{i},2);
end
minerr=inf(1,num);
idx=zeros(p,num);
A=zeros(p,num);
rec=zeros(d,num);

for setcnt=2:2^p-1
    setmask=logical(de2bi(setcnt,p));
    q=sum(setmask);
    % We treat the q=1 case separately at the end
    if q==1
        continue;
    end
    Ni=N(setmask);
    I=zeros(1,q);
    pt=1;
    M=Ni-1;
    setmask=find(setmask);
    E=zeros(d,q);
    while true
        %if q==p
        %    I
        %end
        for i=1:q
            E(:,i)=L{setmask(i)}(:,I(i)+1);
        end
        [y,a]=plane_project2(x,E);  
        mp=find(~(sum(a<0)>0));
        for i=1:numel(mp)
            err=norm(y(:,mp(i))-x(:,mp(i)));
            if err<minerr(mp(i))
                minerr(mp(i))=err;
                idx(:,mp(i))=zeros(p,1);
                idx(setmask,mp(i))=I+1;
                A(:,mp(i))=zeros(p,1);
                A(setmask,mp(i))=a(:,mp(i));
                rec(:,mp(i))=y(:,mp(i));
            end
        end
            
        % generate next set indices
        pt=find(I<M);
        if isempty(pt)
            break;
        end
        pt=pt(1);
        I(pt)=I(pt)+1;
        I(1:pt-1)=0;
    end
end

% q=1 case
for i=1:p
    for j=1:N(i)
        for k=1:num
            err=norm(x(:,k)-L{i}(:,j));
            if err<minerr(k)
                minerr(k)=err;
                idx(:,k)=zeros(p,1);
                idx(i,k)=j;
                A(:,k)=zeros(p,1);
                A(i,k)=1;
                rec(:,k)=L{i}(:,j);
            end
        end
    end
end

end


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
