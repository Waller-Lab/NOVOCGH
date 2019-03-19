close all;clear all;clc;    %Initialize matlab

%This code performs NOVO_HOLO with dual binary masks (0,1) for intensity
%and penalty



% Specify system parameters here
System.verbose=1;           % 1 or 0    Set this value to 1 to display activity, 0 otherwise
System.lambda = 1.030e-6;   % meters    Wavelength of the light
System.focal_SLM = 0.2;     % meters    focal length of the telescope lens after slm.
System.psSLM = 20e-6;       % meters    SLM pixel dimensions
System.Nx = 300;            % int       Number of pixels in X direction
System.Ny = 300;            % int       Number of pixels in Y direction
System.useGPU = 1;          % 1 or 0    Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800). When Nx, Ny is 300x300, no difference from CPU
System.maxiter = 100;        % int       Number of iterations (for all methods explored)
System.GSoffset = 0.01;     % float>0   Regularization constant to allow low light background in 3D Gerchberg Saxton algorithms[ Images, Depths ] = function_dummydata( System );

[ Images, Depths ] = function_dummydata( System );%MODIFY BY YI

KickImages = zeros(size(Images));
KickImages(:,:,2)=Images(:,:,5);
KickImages(:,:,4)=Images(:,:,1);

%Set factor ten proportionally to how much you want to avoid light in a
%specific area... 1 being equally important as having light in other desired spots.  Here, we want no light within the F letter.
priority = 10; %e.g. 0.1 barely eliminates light, 1 somehow, and 10, really kicks light away)
% KickImages(:,:,6) = priority*Images(:,:,6); 
% % and of course erase that letter from the intensity request
% Images(:,:,6) = 0;   

% Relax constraints on letter D. We don't want D, but we don't mind if
% light passes through this area.
Images(:,:,4) = 0;
Depths = Depths/2;


%% Upgrade data format for GPU use
if System.useGPU ==1
    Masks = gpuArray(Images); KickMasks = gpuArray(KickImages);
else    
    Masks = Images; KickMasks = KickImages;    
end
%We also define a coherent light source with unitary intensity .
System.intensity = 1;
System.source = sqrt(System.intensity)*(1/(System.Nx* System.Ny))*ones(System.Nx, System.Ny);

%Make stack of fresnel propagators corresponding to all considered depth levels.
[ HStacks ] = function_Hstacks( System,Depths );


%Display masks and cost function 
% f = figure(1); 
% for i = 1:5
% subplot(2,5,i);    imagesc(squeeze(Images(:,:,i))); axis off; axis image; caxis([0 2]);
% subplot(2,5,5+i);    imagesc(squeeze(KickImages(:,:,i))); axis off; axis image; caxis([0 2]);    
% end
% title('Top row is requested intensity, bottom is penalty for inadvertent illumination')
 


%% Generate hologram by superposition method
[ Superposition ] = function_Superposition( System,HStacks,Masks );


%% Generate hologram by the sequential Gerchberg Saxtion method.
% [ GS ] = function_sequentialGS(System, HStacks, Masks );
[GS]=function_globalGS(System, HStacks, Masks );

% Set up threshold values for both low-pass and high-pass norms, and then compile binary 3D hologram with
% NOVO_CGH algorithm.
NormOptions.HighThreshold = 0.5;
NormOptions.LowThreshold = 0.1;
[NOVOCGH,history] = function_NOVO_CGH_dualbinary( System, HStacks, Masks,KickMasks,Depths,NormOptions );

%We compute the intensity pattern
[  Superposition.IntensityStack] = function_VolumeIntensity( System, Superposition.phase,HStacks );
[  GS.IntensityStack] = function_VolumeIntensity( System, GS.phase,HStacks );
[  NOVOCGH.IntensityStack] = function_VolumeIntensity( System, NOVOCGH.phase,HStacks );

%Display reconstructed intensity at each mask level
% g = figure(2);for i = 1:numel(Depths);
%     imagesc(squeeze(NOVOCGH.IntensityStack(:,:,i)));
%     colormap gray;
%     title(sprintf('Distance z %d', Depths(i)));
%     caxis([0, max(NOVOCGH.IntensityStack(:))]);
%     pause(10/numel(Depths));
% end;
% close(g);



GS.Scores = function_ScoresI(KickMasks, Masks, GS.IntensityStack);
Superposition.Scores = function_ScoresI(KickMasks, Masks, Superposition.IntensityStack);
NOVOCGH.Scores = function_ScoresI(KickMasks, Masks, NOVOCGH.IntensityStack);


% g = figure(2);
% subplot(1,3,1);
% plot(100*Superposition.Scores.PositiveError); hold on;
% plot(100*GS.Scores.PositiveError); hold on;
% plot(100*NOVOCGH.Scores.PositiveError); 
% legend({'Superposition','GS','NOVOHOLO'}); title('Amount of light going through bright area %');
% xlabel('Z depth'); ylabel('Fraction of total intensity (%)');
% subplot(1,3,2);
% plot(100*Superposition.Scores.NegativeError); hold on;
% plot(100*GS.Scores.NegativeError); hold on;
% plot(100*NOVOCGH.Scores.NegativeError); 
% legend({'Superposition','GS','NOVOHOLO'}); title('Amount of light going through dark area %');
% xlabel('Z depth'); ylabel('Fraction of total intensity (%)');
% subplot(1,3,3);
% plot(100*Superposition.Scores.TargetMatcherror); hold on;
% plot(100*GS.Scores.TargetMatcherror); hold on;
% plot(100*NOVOCGH.Scores.TargetMatcherror); 
% legend({'Superposition','GS','NOVOHOLO'}); title('Error in Intensity distribution');
% xlabel('Z depth'); ylabel('Error (%)');

% save('Workspace_NOVO_CGH_DualBinary.mat');
% figure();plot(history.fval);
% xlim([1 (System.maxiter+1)]);
% grid on;
% set(gca,'fontsize',16);
% xlabel('iteration');ylabel('cost function');
% saveas(gcf,'history_costfunction.fig','fig');
