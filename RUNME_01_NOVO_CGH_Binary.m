close all;clear all;clc;    %Initialize matlab

%This code performs NOVO_HOLO with simple binary masks (1 means no light, 0 means avoid light)

% Specify system parameters here
System.verbose=1;           % 1 or 0    Set this value to 1 to display activity, 0 otherwise
System.lambda = 1.030e-6;   % meters    Wavelength of the light
System.focal_SLM = 0.2;     % meters    focal length of the telescope lens after slm.
System.psSLM = 20e-6;       % meters    SLM pixel dimensions
System.Nx = 300;            % int       Number of pixels in X direction
System.Ny = 300;            % int       Number of pixels in Y direction
System.useGPU = 1;          % 1 or 0    Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800).
System.maxiter = 100;        % int       Number of iterations (for all methods explored)
System.GSoffset = 0.01;     % float>0   Regularization constant to allow low light background in 3D Gerchberg Saxton algorithms

%First we generate target hologram,s Masks, and associated depths, Depths,
%with aproporiate format (if GPU used) from a dummy set of images, Images.
% using our dummy data generator
[ Images, Depths ] = function_dummydata( System );
%Images = Images(:,:,[2 4 6]); Depths = Depths([2 4 6]);

Depths = Depths/2;
if System.useGPU ==1
    Masks = gpuArray(Images);
else
    Masks = Images;
end;


%We also define a coherent light source with unitary intensity .
System.intensity = 1;
System.source = sqrt(System.intensity)*(1/(System.Nx* System.Ny))*ones(System.Nx, System.Ny);

%Make stack of fresnel propagators corresponding to all considered depth levels.
[ HStacks ] = function_Hstacks( System,Depths );

%% Generate hologram by superposition method
[ Superposition ] = function_Superoposition( System,HStacks,Masks );


%% Generate hologram by the sequential Gerchberg Saxtion method.
[ GS ] = function_sequentialGS(System, HStacks, Masks );


% Set up threshold values for both low-pass and high-pass norms, and then compile binary 3D hologram with
% NOVO_CGH algorithm.
NormOptions.HighThreshold = 0.5;
NormOptions.LowThreshold = 0.1;
[ NOVOCGH.hologram,NOVOCGH.phase ] = function_NOVO_CGH_binary( System, HStacks, Masks,Depths,NormOptions );

%We compute the intensity pattern
[  Superposition.IntensityStack] = function_VolumeIntensity( System, Superposition.phase,HStacks );
[  GS.IntensityStack] = function_VolumeIntensity( System, GS.phase,HStacks );
[  NOVOCGH.IntensityStack] = function_VolumeIntensity( System, NOVOCGH.phase,HStacks );

%Display reconstructed intensity at each mask level
f = figure(1);for i = 1:numel(Depths);
    imagesc(squeeze(NOVOCGH.IntensityStack(:,:,i)));
    colormap gray;
    title(sprintf('Distance z %d', Depths(i)));
    caxis([0, max(NOVOCGH.IntensityStack(:))]);
    pause(5/numel(Depths));
end;
close(f);



GS.Scores = function_ScoresI(1-Masks, Masks, GS.IntensityStack);
Superposition.Scores = function_ScoresI(1-Masks, Masks, Superposition.IntensityStack);
NOVOCGH.Scores = function_ScoresI(1-Masks, Masks, NOVOCGH.IntensityStack);


g = figure(2);
subplot(1,3,1);
plot(100*Superposition.Scores.PositiveError); hold on;
plot(100*GS.Scores.PositiveError); hold on;
plot(100*NOVOCGH.Scores.PositiveError); 
legend({'Superposition','GS','NOVOHOLO'}); title('Amount of light going through bright area %');
xlabel('Z depth'); ylabel('Fraction of total intensity (%)');
subplot(1,3,2);
plot(100*Superposition.Scores.NegativeError); hold on;
plot(100*GS.Scores.NegativeError); hold on;
plot(100*NOVOCGH.Scores.NegativeError); 
legend({'Superposition','GS','NOVOHOLO'}); title('Amount of light going through dark area %');
xlabel('Z depth'); ylabel('Fraction of total intensity (%)');
subplot(1,3,3);
plot(100*Superposition.Scores.TargetMatcherror); hold on;
plot(100*GS.Scores.TargetMatcherror); hold on;
plot(100*NOVOCGH.Scores.TargetMatcherror); 
legend({'Superposition','GS','NOVOHOLO'}); title('Intensity_error');
xlabel('Z depth'); ylabel('Error (%)');

