close all;
clear all;

tag = '1024by1272_ABCDEFGHIJ';


%% Setup params
% All length value has unit meter in this file.
% The 3d region is behind lens after SLM.

resolutionScale = 1; % The demagnification scale of tubelens and objective. f_tube/f_objective
lambda = 0.633e-6;  % Wavelength
focal_SLM = 0.2; % focal length of the lens after slm.
psSLM = 12.5e-6;      % Pixel Size (resolution) at the scattered 3D region
Nx = 1272;       % Number of pixels in X direction
Ny = 1024;       % Number of pixels in Y direction

psXHolograph = lambda * focal_SLM/ psSLM / resolutionScale / Nx;      % Pixel Size (resolution) at the scattered 3D region
psYHolograph = lambda * focal_SLM/ psSLM / resolutionScale / Ny;      % Pixel Size (resolution) at the scattered 3D region

useGPU = 1;     % Use GPU to accelerate computation. Effective when Nx, Ny is large (e.g. 600*800).
z = [-100 : 2: 100] * 16e-4 ;   % Depth level requested in 3D region.




%% Complex Target

load('ABCDEFGHIJK');
n = 10;
% % Preprocess image
resizeX =450;
resizeY = 600;
cropxmin = 1; cropxmax=450; cropymin=1; cropymax=600;
shiftx = -50;
shifty = 0;
padx = 80;
pady = 80;
masks = imresize(masks, [resizeX, resizeY]);
masks = circshift(masks(cropxmin:cropxmax, cropymin:cropymax, :), [shiftx, shifty, 0]);
masks = padarray(masks, [padx, pady, 0]);
ims = imresize(permute(masks(:,:,1:n), [2,1,3]) ,[Nx, Ny]);


step = floor(50 / n) + 1;
indices = [21:step:71];
if numel(indices) < n
    indices(end+1) = indices(end) + step;
end

imdepths = z(indices);
z = z(indices);
nfocus = floor(numel(z)/2)+1;
focus = (z(end) + z(1))/2;

zdisplay = imdepths;

%% Setup binary target intensity matrix named 'Masks'

Masks = zeros(Nx, Ny, numel(z));
if useGPU
    Masks = gpuArray(Masks);
end

figure()
for i = 1 : numel(z)
    Masks(:,:,i) = generateComplexMask( z(i), Nx, Ny, ims, imdepths);
    imagesc(Masks(:,:,i));pause(5/numel(z))
end
Masks(Masks >= 0.5) = 1;
Masks(Masks < 0.5) = 0;
maskfun = @(i1, i2) Masks(:,:,i1:i2);



%% Precompute Fresnel Kernel Function
HStacks = zeros(Nx, Ny, numel(z));

for i = 1 : numel(z)
    HStacks(:,:,i) = GenerateFresnelPropagationStack(Nx, Ny, z(i)-focus, lambda, psXHolograph,psYHolograph, 0);
end
if useGPU
    kernelfun = @(i1, i2) gpuArray(HStacks(:,:,i1:i2));
else
    kernelfun = @(i1, i2) HStacks(:,:,i1:i2);
end


%% Pick Source Initialization method

%This sets a coherent light source.
intensity = 1;
source = sqrt(intensity) * ones(Nx, Ny);
tag = [tag, '_uniform'];
source = source/(max(max(source)));
tic;

%% Superposition


tic
hologram = zeros(Nx, Ny);
for i = 1:numel(imdepths)
    HStack = GenerateFresnelPropagationStack(Nx, Ny, imdepths(i)-focus, lambda, psXHolograph,psYHolograph, useGPU);
    target = Masks(:,:,i);
    target_phase = randn(Nx, Ny);
    imagez = sqrt(target) .* exp(1i * target_phase);
    im =  ifft2(ifftshift(imagez))./HStack;
    hologram = hologram + im;
end
toc

phase = gather(angle(hologram));

%% plot
Ividmeas = zeros(Nx, Ny, numel(zdisplay));
usenoGPU = 0;
figure();
phase = mod(phase, 2*pi) - pi;
for i = 1:numel(zdisplay)
    HStack = GenerateFresnelPropagationStack(Nx,Ny,zdisplay(i) -focus, lambda, psXHolograph,psYHolograph, usenoGPU);
    imagez = fresnelProp(phase, source, HStack);
    Ividmeas(:,:,i) = imagez;
    imagesc(imagez);colormap gray;title(sprintf('Distance z %d', zdisplay(i)));
    maxInten = sum(sum(source.^2))*500;
    colorbar;
    caxis([0, 5e8]);
    pause(5/numel(zdisplay));
end
pct9999 = prctile(Ividmeas(:), 99.99);
prctile(Ividmeas(:), 99.9);
pct100 = max(Ividmeas(:));


phase = mod(phase, 2*pi) - pi;
x0 = phase;
save(['superposition' tag '.mat'], 'phase', 'source', 'z');





%% sequential GS


ims = Masks;
maxiter = 100;
figure();
imageInten = gpuArray(zeros(Nx, Ny, numel(imdepths)));

tic;

mu = 0.01;
im = source.*exp(1i * x0);

for n = 1:maxiter+1
    if mod(n, 100) == 0
        display(n);
    end
    index = 1:numel(imdepths);
    tempim = 0 * im;
    for i = index
        HStack = GenerateFresnelPropagationStack(Nx, Ny, imdepths(i)-focus, lambda, psXHolograph,psYHolograph, useGPU);
        imagez = fftshift(fft2(im .* HStack));
        target = ims(:,:,i) + mu;
        imagez = sqrt(target) .* exp(1i * angle(imagez));
        tempim = tempim + ifft2(ifftshift(imagez))./HStack;
    end
    
    im = source.*exp(1i * angle(tempim));
    
end

toc;

phase = gather(angle(im));
%% plot
Ividmeas = zeros(Nx, Ny, numel(zdisplay));
usenoGPU = 0;
figure();
phase = mod(phase, 2*pi) - pi;
for i = 1:numel(zdisplay)
    HStack = GenerateFresnelPropagationStack(Nx,Ny,zdisplay(i) - focus, lambda, psXHolograph,psYHolograph, usenoGPU);
    imagez = fresnelProp(phase, source, HStack);
    Ividmeas(:,:,i) = imagez;
    imagesc(imagez);colormap gray;title(sprintf('Distance z %d', zdisplay(i)));
    maxInten = sum(sum(source.^2))*500;
    colorbar;
    caxis([0, 5e8]);
    pause(5/numel(zdisplay));
end
phase = mod(phase, 2*pi) - pi;
save(['Sequential_gs' tag '.mat'], 'phase', 'source', 'z');


%% Optimization

thresholdh = 0.5 * (pct9999 + pct100)/2;
thresholdl = 0.1 * thresholdh;
maxiter = 50;

f = @(x)FunObj( x, source, z, Nx, Ny, thresholdh, thresholdl, maskfun, kernelfun, useGPU);
tic;
matlab_options = optimoptions('fmincon','GradObj','on', 'display', 'iter', ...
    'algorithm','interior-point','Hessian','lbfgs', 'MaxFunEvals', 500, 'MaxIter', maxiter,...
    'TolX', 1e-20, 'TolFun', 1e-12);
lb = -inf(Nx*Ny, 1);
ub = inf(Nx*Ny, 1);
nonlcon = [];
phase = fmincon(f,phase,[],[],[],[],lb,ub,nonlcon,matlab_options);
phase = reshape(phase, [Nx, Ny]);
phase = mod(phase, 2*pi) - pi;
hologram = floor(mod(phase, 2*pi)/2/pi * 255);

toc;

%% plot
Ividmeas = zeros(Nx, Ny, numel(zdisplay));
usenoGPU = 0;
figure();
for i = 1:numel(zdisplay)
    HStack = GenerateFresnelPropagationStack(Nx,Ny,zdisplay(i) - focus, lambda, psXHolograph,psYHolograph, usenoGPU);
    imagez = fresnelProp(phase, source, HStack);
    
    Ividmeas(:,:,i) = imagez;
    imagesc(imagez);colormap gray;title(sprintf('Distance z %d', zdisplay(i)));colorbar;
    caxis([0, 3e8]);
    pause(5/numel(zdisplay));
end

prctile(Ividmeas(:), 99.99)
save(['optimization_' tag '.mat'], 'source', 'phase', ...
    'thresholdl', 'thresholdh', 'z', 'resizeX', 'resizeY', 'shiftx', 'shifty');
