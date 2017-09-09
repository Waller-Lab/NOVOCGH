close all;clear all;clc;    %Initialize matlab

% Specify system parameters here
System.verbose=1;           % 1 or 0    Set this value to 1 to display activity, 0 otherwise
System.lambda = 1.030e-6;   % meters    Wavelength of the light
System.focal_SLM = 0.2;     % meters    focal length of the telescope lens after slm.
System.psSLM = 20e-6;       % meters    SLM pixel dimensions
System.Nx = 300;            % int       Number of pixels in X direction
System.Ny = 300;            % int       Number of pixels in Y direction
System.useGPU = 1;          % 1 or 0    Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800).
System.maxiter = 50;        % int       Number of iterations (for all methods explored)
System.GSoffset = 0.01;     % float>0   Regularization constant to allow low light background in 3D Gerchberg Saxton algorithms

% Specify Low and High threshold for threshold-based cost functions
NormOptions.HighThreshold = 0.5;
NormOptions.LowThreshold = 0.1;

%Specify Illumination pattern at the SLM, Uniform here, but tunable in
%general.
System.intensity = 1;
System.source = sqrt(System.intensity)*(1/(System.Nx* System.Ny))*ones(System.Nx, System.Ny);

% Generate Dummy Data, and format for GPU use if needed.
[ TargetIntensity, TargetDarkenss , Depths ] = function_dummydata( System );
[LX, LY, LZ] = size(TargetIntensity); 
if System.useGPU ==1
    Masks = gpuArray(TargetIntensity);
    KickMasks = gpuArray(TargetDarkenss);
else
    Masks = TargetIntensity;
    KickMasks = TargetDarkenss;
end;

%Compute the kernels for propagation at all defined depths
[ HStacks ] = function_Hstacks( System,Depths );

%Compute holograms with Iterative methods
Hologram{1} = function_Superposition(System, HStacks, Masks );
Hologram{2} = function_sequentialGS(System, HStacks, Masks );
Hologram{3} = function_globalGS(System, HStacks, Masks );
%Compute holograms with various implementations of NOVO-CGH
% 1 - Binary targets 
Hologram{4} = function_NOVO_CGH_binary( System, HStacks, Masks,Depths,NormOptions );

%2 - Bright voxels, Dark voxels, and relaxed constraints anywhere where KickMasks = Masks = 0, 
Hologram{5} = function_NOVO_CGH_dualbinary( System, HStacks, Masks, KickMasks,Depths,NormOptions );

%3 - Variable Intensity Holograms with a threshold based cost function
Hologram{6} = function_NOVO_CGH_VarI( System, HStacks, Masks,Depths,NormOptions);

%4 - Variable Intensity hologram with a Euclidean cost function
Hologram{7} = function_NOVO_CGH_VarIEuclid( System, HStacks, Masks,Depths);

%5 - Two-photon hologram (specify the desired two photon absorption pattern directly 
Hologram{8} = function_NOVO_CGH_TPEuclid( System, HStacks, Masks.^2,Depths);


% Compute intensity at each specified z plane for all holograms
for i = 1:7
Display{i} = function_VolumeIntensity( System, Hologram{i}.phase,HStacks );
end
Display{8} = function_VolumeIntensity( System, Hologram{i}.phase,HStacks ).^2; % For two photon repsonse
Display{9} = double(Masks>0);
Display{10} = double(KickMasks>0);


% Here we just rescale the global intensity in all volume images to best match the target
% intensity in the voxels where V+>0
for i = 1:8
    pick = Display{i};
    pick =  pick*(sum(pick(:).*gather(Masks(:)))/sum(pick(:).^2));   
    Display{i} = pick;
end

% We display first the holograms where only V+ is definesd then the V+
% target , the dark voxels (V-)
% and finally the holograms where both V+ and V- are used together for 3 level
% control (Bright, Dark, and don't care).
OrderPlots = [1 2 3 4 6 7 8 9 10 5];
Holograms_typesLine1 = {'Superposition', 'Sequential', 'Global', 'NOVO-CGH', 'NOVO-CGH', 'NOVO-CGH', 'NOVO-CGH', 'NOVO-CGH', 'Target', 'Target '};
Holograms_typesLine2 = {'Superposition', 'GS', 'GS', 'NOVO-Binary', 'Binary+dark ROIs', 'Variable I', 'Euclidian norm', 'Two-photon CGH', 'Intensity', 'Dark Voxels'};

%Making the final figure
f = figure(1)
for i = 1:10
    for j = 1:LZ
        subplot(LZ,10,i+(j-1)*10)       
        imagesc(gather(squeeze(Display{OrderPlots(i)}(:,:,j))));
        axis off; axis image; colormap hot;
        title({Holograms_typesLine1{OrderPlots(i)}, Holograms_typesLine2{OrderPlots(i)} });
        caxis([0 1.5])
    end
       
end






