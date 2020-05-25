prompt = 'Running this code will clear your workspace. Would you like to continue?[Y/N]';
str = input(prompt,'s');
if (str=='N') || (str=='n')
    ;
elseif (str=='Y') || (str=='y')
    
clear;close all;clc
load('demo.mat')

[Parameters] = BCMParameters(endmembers);

%% 	BCM Unmixing
%%% Approximate running time: 70 seconds per BCM approach*. 
%%% *Measured on a desktop PC with Intel i7 3.20 GHz processor and 12 GB RAM.

disp('BCM-Spectral-QP Unmixing...');
[P1] = BCM(Xim, Parameters, 1);%BCM-Spectral-QP
P1r = reshape(P1,[13 19 4]);

disp('BCM-Spectral-MH Unmixing...');
[P2] = BCM(Xim, Parameters, 2);%BCM-Spectral-MH
P2r = reshape(P2,[13 19 4]);

disp('BCM-Spatial-QP Unmixing...');
[P3] = BCM(Xim, Parameters, 3);%BCM-Spatial-QP
P3r = reshape(P3,[13 19 4]);

disp('BCM-Spatial-MH Unmixing...');
[P4] = BCM(Xim, Parameters, 4);%BCM-Spatial-MH
P4r = reshape(P4,[13 19 4]);
%% FCLS and NCM Unmixing (for comparison)
X = reshape(Xim,[size(Xim,1)*size(Xim,2),size(Xim,3)]);

disp('FCLS Unmixing...');
[ PFCLS ] = hyperFcls( X', Parameters.Emean' );  %FCLS
PFCLS = PFCLS';
PFCLSr = reshape(PFCLS,[13 19 4]);

disp('NCM-QP Unmixing...');
[PestimateNCM1] = unmix2(X', Parameters.Emean'); %NCM-QP
PNCM1r = reshape(PestimateNCM1,[13 19 4]);

disp('NCM-MH Unmixing...');
[PestimateNCM2] = unmixGaussian(X,Parameters);   %NCM-MH
PNCM2r = reshape(PestimateNCM2,[13 19 4]);
%% Plot proportion maps
fig = figure(100);
set(gcf,'Position',get(0,'ScreenSize'));
subplot(2,2,1);imagesc(P1r(:,:,1),[0 1]);title('Asphalt');axis equal tight
subplot(2,2,2);imagesc(P1r(:,:,2),[0 1]);title('Yellow Curb');axis equal tight
subplot(2,2,3);imagesc(P1r(:,:,3),[0 1]);title('Grass');axis equal tight
subplot(2,2,4);imagesc(P1r(:,:,4),[0 1]);title('Oak Leaves');axis equal tight
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf BCM-Spectral-QP Unmixing: Proportion Maps','HorizontalAlignment' ,'center','VerticalAlignment', 'top')


fig = figure(200);
set(gcf,'Position',get(0,'ScreenSize'));
subplot(2,2,1);imagesc(P2r(:,:,1),[0 1]);title('Asphalt');axis equal tight
subplot(2,2,2);imagesc(P2r(:,:,2),[0 1]);title('Yellow Curb');axis equal tight
subplot(2,2,3);imagesc(P2r(:,:,3),[0 1]);title('Grass');axis equal tight
subplot(2,2,4);imagesc(P2r(:,:,4),[0 1]);title('Oak Leaves');axis equal tight
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf BCM-Spectral-MH Unmixing: Proportion Maps','HorizontalAlignment' ,'center','VerticalAlignment', 'top')


fig = figure(300);
set(gcf,'Position',get(0,'ScreenSize'));
subplot(2,2,1);imagesc(P3r(:,:,1),[0 1]);title('Asphalt');axis equal tight
subplot(2,2,2);imagesc(P3r(:,:,2),[0 1]);title('Yellow Curb');axis equal tight
subplot(2,2,3);imagesc(P3r(:,:,3),[0 1]);title('Grass');axis equal tight
subplot(2,2,4);imagesc(P3r(:,:,4),[0 1]);title('Oak Leaves');axis equal tight
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf BCM-Spatial-QP Unmixing: Proportion Maps','HorizontalAlignment' ,'center','VerticalAlignment', 'top')


fig = figure(400);
set(gcf,'Position',get(0,'ScreenSize'));
subplot(2,2,1);imagesc(P4r(:,:,1),[0 1]);title('Asphalt');axis equal tight
subplot(2,2,2);imagesc(P4r(:,:,2),[0 1]);title('Yellow Curb');axis equal tight
subplot(2,2,3);imagesc(P4r(:,:,3),[0 1]);title('Grass');axis equal tight
subplot(2,2,4);imagesc(P4r(:,:,4),[0 1]);title('Oak Leaves');axis equal tight
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf BCM-Spatial-MH Unmixing: Proportion Maps','HorizontalAlignment' ,'center','VerticalAlignment', 'top')


fig = figure(500);
set(gcf,'Position',get(0,'ScreenSize'));
subplot(2,2,1);imagesc(PFCLSr(:,:,1),[0 1]);title('Asphalt');axis equal tight
subplot(2,2,2);imagesc(PFCLSr(:,:,2),[0 1]);title('Yellow Curb');axis equal tight
subplot(2,2,3);imagesc(PFCLSr(:,:,3),[0 1]);title('Grass');axis equal tight
subplot(2,2,4);imagesc(PFCLSr(:,:,4),[0 1]);title('Oak Leaves');axis equal tight
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf FCLS Unmixing: Proportion Maps','HorizontalAlignment' ,'center','VerticalAlignment', 'top')


fig = figure(600);
set(gcf,'Position',get(0,'ScreenSize'));
subplot(2,2,1);imagesc(PNCM1r(:,:,1),[0 1]);title('Asphalt');axis equal tight
subplot(2,2,2);imagesc(PNCM1r(:,:,2),[0 1]);title('Yellow Curb');axis equal tight
subplot(2,2,3);imagesc(PNCM1r(:,:,3),[0 1]);title('Grass');axis equal tight
subplot(2,2,4);imagesc(PNCM1r(:,:,4),[0 1]);title('Oak Leaves');axis equal tight
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf NCM-QP Unmixing: Proportion Maps','HorizontalAlignment' ,'center','VerticalAlignment', 'top')


fig = figure(700);
set(gcf,'Position',get(0,'ScreenSize'));
subplot(2,2,1);imagesc(PNCM2r(:,:,1),[0 1]);title('Asphalt');axis equal tight
subplot(2,2,2);imagesc(PNCM2r(:,:,2),[0 1]);title('Yellow Curb');axis equal tight
subplot(2,2,3);imagesc(PNCM2r(:,:,3),[0 1]);title('Grass');axis equal tight
subplot(2,2,4);imagesc(PNCM2r(:,:,4),[0 1]);title('Oak Leaves');axis equal tight
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf NCM-MH Unmixing: Proportion Maps','HorizontalAlignment' ,'center','VerticalAlignment', 'top')

%% Compute PError against manual ground truth and display results
load('Ptrue.mat')
N = size(Xim,1)*size(Xim,2);
[~, PErrorOutputmean1, ~] = PError(N,Ptrue,P1);
[~, PErrorOutputmean2, ~] = PError(N,Ptrue,P2);
[~, PErrorOutputmean3, ~] = PError(N,Ptrue,P3);
[~, PErrorOutputmean4, ~] = PError(N,Ptrue,P4);
[~, PErrorOutputmeanF, ~] = PError(N,Ptrue,PFCLS);
[~, PErrorOutputmeanN1, ~] = PError(size(X,1),Ptrue,PestimateNCM1);
[~, PErrorOutputmeanN2, ~] = PError(size(X,1),Ptrue,PestimateNCM2);


% Display PErrorOutputmean result
fprintf('PErrorOutputMean for FCLS: ');
fprintf(' %f \n', PErrorOutputmeanF);

fprintf('PErrorOutputMean for NCM-QP: ');
fprintf(' %f \n', PErrorOutputmeanN1);

fprintf('PErrorOutputMean for NCM-MH: ');
fprintf(' %f \n', PErrorOutputmeanN2);

fprintf('PErrorOutputMean for BCM-Spectral-QP: ');
fprintf(' %f \n', PErrorOutputmean1);

fprintf('PErrorOutputMean for BCM-Spectral-MH: ');
fprintf(' %f \n', PErrorOutputmean2);

fprintf('PErrorOutputMean for BCM-Spatial-QP: ');
fprintf(' %f \n', PErrorOutputmean3);

fprintf('PErrorOutputMean for BCM-Spatial-MH: ');
fprintf(' %f \n', PErrorOutputmean4);

end
