function [index_set, abundances, reconstruct, err] = AAM_no_subset( X, L, K )
% AAM_NO_SUBSET Alternating angle minimization implementation of the MESMA
% algorithm, called by AAM
%
% This function executes the AAM algorithm on a set of libraries. For every
% pixel in x, an endmember set is constructed containing a single endmember
% from each library. The indices and abundances are returned. The algorithm
% employed is iterative alternating angle minimization.
%
% Input:  X: the N data points to unmix in d dimensions, (d,N)
%         L: spectral library, cell array of p elements of size (d,~)
%         K: Optional number of iterations. Default 3
% Output: index_set identifying the endmembers from each library, (p,N)
%         abundances with respect to these endmembers, (p,N)
%         reconstruct contains the reconstructed spectral (d,N)
%         err contains the reconstruction error (Euclidean distance) (1,N)
%
%
% Rob Heylen, 2016, University of Antwerp.



% Turn this on to activate sanity checks on the libraries and input values.
% Turn this off if you are sure there are no doubles in the libraries, and
% you do not want to check if pixels are contained in the libraries.
library_check=1;

% Initializations
if nargin==2
    K=3;
end
[d,numpx]=size(X);
p=numel(L);
for i=1:p
    N(i)=size(L{i},2);
end
flag=0;
index_set=zeros(p,numpx);
abundances=zeros(p,numpx);
F=zeros(d,p);
I=ones(p,1);
reconstruct=zeros(d,numpx);
err=zeros(1,numpx);


% Main loop over all pixels
for px=1:numpx
    x=X(:,px);
    
    % Check if x is a library member. If so, we can finish immediately
    if library_check
        for i=1:p
            %if sum(sum(abs(L{i}-x*ones(1,N(i))))==0)>0
            if numel(find(~sum(abs(L{i}-x*ones(1,N(i))))))>0
                I=ones(p,1);
                I(i)=find(sum(abs(L{i}-x*ones(1,N(i))))==0);
                index_set(:,px)=I;
                flag=1;
                break;
            end
        end
    end
    if flag
        flag=0;
        continue;
    end

    % Create random initial endmember set    
    for i=1:p
        I(i)=ceil(rand*N(i));
        F(:,i)=L{i}(:,I(i));
    end
    
    % Iterate K times
    for it=1:K
        % Alternating angle optimization
        for i=1:p
            % Calculate angles
            Fi=F(:,[1:i-1 i+1:p]);       % Pivot
            Gi=[Fi x];                   % Plane through pivot and pixel
            E1=plane_project2(L{i},Fi);  % Project library onto pivot
            E2=plane_project2(L{i},Gi);  % Project library onto plane
            p1=sqrt(sum((E1-L{i}).^2));  % Distances from library to pivot
            p2=sqrt(sum((E2-L{i}).^2));  % Distances from library to plane
            ang=asin(p2./p1);            % Resulting angles
            
            % Find angles that should be inverted
            mask=(x-plane_project2(x,Fi))'*(E2-E1)<0;
            ang(mask)=pi-ang(mask);
            
            % Identify minimal angle, update index set
            [~,I(i)]=min(ang);
            
            % Update endmember set
            F(:,i)=L{i}(:,I(i));
        end
    end
    
    % Update index_set with obtained indices
    index_set(:,px)=I;
end

% Unmixing phase
E=zeros(d,p);
go=0;
for px=1:numpx
    for i=1:p
        E(:,i)=L{i}(:,index_set(i,px));
    end
    % Plug in your favorite unmixing program here
    
    if go==0
        [at,opt]=FCLSU_fast2(X(:,px)',E);
        abundances(:,px)=at';
        go=1;
    else
        at=FCLSU_fast2(X(:,px)',E,opt);
        abundances(:,px)=at';
    end
    
    % Reconstruction
    reconstruct(:,px)=E*abundances(:,px);
    
    % Error
    err(px)=norm(reconstruct(:,px)-X(:,px));
end

end


