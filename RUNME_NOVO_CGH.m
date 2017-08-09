close all;clear all;clc;    %Initialize matlab

% Specify system parameters here
System.verbose=1;           % 1 or 0    Set this value to 1 to display activity, 0 otherwise
System.lambda = 1.030e-6;   % meters    Wavelength of the light
System.focal_SLM = 0.2;     % meters    focal length of the telescope lens after slm.
System.psSLM = 20e-6;       % meters    SLM pixel dimensions
System.Nx = 300;            % int       Number of pixels in X direction
System.Ny = 300;            % int       Number of pixels in Y direction
System.useGPU = 1;          % 1 or 0    Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800).
System.maxiter = 10;        % int       Number of iterations (for all methods explored)
System.GSoffset = 0.01;     % float>0   Regularization constant to allow low light background in 3D Gerchberg Saxton algorithms

%First we generate target hologram,s Masks, and associated depths, Depths,
%with aproporiate format (if GPU used) from a dummy set of images, Images.
% using our dummy data generator
[ Images, Depths ] = function_dummydata( System );
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
[ Superposition.hologram ] = function_Superoposition( System,HStacks,Masks );


%% Generate hologram by the sequential Gerchberg Saxtion method.
[ GS.hologram ] = function_sequentialGS(System, HStacks, Masks );
GS.phase = gather(angle(GS.hologram));

% Set up threshold values for both low-pass and high-pass norms, and then compile binary 3D hologram with
% NOVO_CGH algorithm.
NormOptions.HighThreshold = 0.5;
NormOptions.LowThreshold = 0.1;
[ NOVOCGH.hologram,NOVOCGH.phase ] = function_NOVO_CGH_binary( System, HStacks, Masks,Depths,NormOptions );

%We compute the intensity pattern
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

