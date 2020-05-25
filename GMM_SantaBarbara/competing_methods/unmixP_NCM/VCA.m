function [Ae, indice, Rp] = vca(R,varargin)

% Vertex Component Analysis
%
% [Ae, indice, Rp ]= vca(R,'Endmembers',p,'SNR',r,'verbose',v)
%
% ------- Input variables -------------
%  R - matrix with dimensions L(channels) x N(pixels)
%      each pixel is a linear mixture of p endmembers
%      signatures R = M x s, where s = gamma x alfa
%      gamma is a illumination perturbation factor and
%      alfa are the abundance fractions of each endmember.
% 'Endmembers'
%          p - positive integer number of endmembers in the scene
%
% ------- Output variables -----------
% A - estimated mixing matrix (endmembers signatures)
% indice - pixels that were chosen to be the most pure
% Rp - Data matrix R projected.   
%
% ------- Optional parameters---------
% 'SNR'
%          r - (double) signal to noise ratio (dB)
% 'verbose'
%          v - [{'on'} | 'off']
% ------------------------------------
%
% Authors: José Nascimento (zen@isel.pt) 
%          José Bioucas Dias (bioucas@lx.it.pt) 
% Copyright (c)
% version: 2.1 (7-May-2004)
%
% For any comment contact the authors
%
% more details on:
% José M. P. Nascimento and José M. B. Dias 
% "Vertex Component Analysis: A Fast Algorithm to Unmix Hyperspectral Data"
% submited to IEEE Trans. Geosci. Remote Sensing, vol. .., no. .., pp. .-., 2004
% 
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         verbose = 'on'; % default
         snr_input = 0;  % default this flag is zero,
                         % which means we estimate the SNR
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking for input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         dim_in_par = length(varargin);
         if (nargin - dim_in_par)~=1
            error('Wrong parameters');
         elseif rem(dim_in_par,2) == 1
            error('Optional parameters should always go by pairs');
         else
            for i = 1 : 2 : (dim_in_par-1)
                switch lower(varargin{i})
                  case 'verbose'
                       verbose = varargin{i+1};
                  case 'endmembers'     
                       p = varargin{i+1};
                  case 'snr'     
                       SNR = varargin{i+1};
                       snr_input = 1;       % flag meaning that user gives SNR 
                  otherwise
                       fprintf(1,'Unrecognized parameter:%s\n', varargin{i});
                end %switch
            end %for
         end %if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         if isempty(R)
            error('there is no data');
         else
            [L N]=size(R);  % L number of bands (channels)
                            % N number of pixels (LxC) 
         end                   
               
         if (p<0 | p>L | rem(p,1)~=0),  
            error('ENDMEMBER parameter must be integer between 1 and L');
         end
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         if snr_input==0,
            r_m = mean(R,2);      
            R_m = repmat(r_m,[1 N]); % mean of each band
            R_o = R - R_m;           % data with zero-mean 
            [Ud,Sd,Vd] = svds(R_o*R_o'/N,p);  % computes the p-projection matrix 
            x_p =  Ud' * R_o;                 % project the zero-mean data onto p-subspace
            
            SNR = estimate_snr(R,r_m,x_p);
            
            if strcmp (verbose, 'on'), fprintf(1,'SNR estimated = %g[dB]\n',SNR); end
         else   
            if strcmp (verbose, 'on'), fprintf(1,'input    SNR = %g[dB]\t',SNR); end
         end

         SNR_th = 15 + 10*log10(p);
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choosing Projective Projection or 
%          projection to p-1 subspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         if SNR < SNR_th,   
                if strcmp (verbose, 'on'), fprintf(1,'... Select the projective proj.\n',SNR); end
                
                d = p-1;
                if snr_input==0, % it means that the projection is already computed
                     Ud= Ud(:,1:d);    
                else
                     r_m = mean(R,2);      
                     R_m = repmat(r_m,[1 N]); % mean of each band
                     R_o = R - R_m;           % data with zero-mean 
         
                     [Ud,Sd,Vd] = svds(R_o*R_o'/N,d);  % computes the p-projection matrix 

                     x_p =  Ud' * R_o;                 % project thezeros mean data onto p-subspace

                end
                
                Rp =  Ud * x_p(1:d,:) + repmat(r_m,[1 N]);      % again in dimension L
                
                x = x_p(1:d,:);             %  x_p =  Ud' * R_o; is on a p-dim subspace
                c = max(sum(x.^2,1))^0.5;
                y = [x ; c*ones(1,N)] ;
         else
                if strcmp (verbose, 'on'), fprintf(1,'... Select proj. to p-1\n',SNR); end
             
                d = p;
                [Ud,Sd,Vd] = svds(R*R'/N,d);         % computes the p-projection matrix 
                
                x_p = Ud'*R;
                Rp =  Ud * x_p(1:d,:);      % again in dimension L (note that x_p has no null mean)
                
                x =  Ud' * R;
                u = mean(x,2);        %equivalent to  u = Ud' * r_m
                y =  x./ repmat( sum( x .* repmat(u,[1 N]) ) ,[d 1]);

          end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VCA algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indice = zeros(1,p);
A = zeros(p,p);
A(p,1) = 1;

for i=1:p
      w = rand(p,1);   
      f = w - A*pinv(A)*w;
      f = f / sqrt(sum(f.^2));
      
      v = f'*y;
      [v_max indice(i)] = max(abs(v));
      A(:,i) = y(:,indice(i));        % same as x(:,indice(i))
end
Ae = Rp(:,indice);

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of the vca function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function snr_est = estimate_snr(R,r_m,x)

         [L N]=size(R);           % L number of bands (channels)
                                  % N number of pixels (Lines x Columns) 
         [p N]=size(x);           % p number of endmembers (reduced dimension)

         P_y = sum(R(:).^2)/N;
         P_x = sum(x(:).^2)/N + r_m'*r_m;
         snr_est = 10*log10( (P_x - p/L*P_y)/(P_y- P_x) );
return;
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
