function [ HStacks ] = function_Hstacks( System,z )
%Precompute Fresnel propagation kernels at depths specified by vector z 
%for fast computation of wave propagation in all algorithms.

if System.verbose == 1
    tic;
    fprintf('Kernel computation begins...');
end;
psXHolograph = System.lambda * System.focal_SLM/ System.psSLM  / System.Nx;      % Pixel Size (resolution) at the scattered 3D region
psYHolograph = System.lambda * System.focal_SLM/ System.psSLM  / System.Ny;      % Pixel Size (resolution) at the scattered 3D region
if System.useGPU ==1
    HStacks = zeros(System.Nx, System.Ny, numel(z), 'gpuArray');
else
    HStacks = zeros(System.Nx, System.Ny, numel(z));
end
for i = 1 : numel(z)
    HStacks(:,:,i) = function_GenerateFresnelPropagationStack(System.Nx, System.Ny, z(i), System.lambda, psXHolograph,psYHolograph, System.useGPU);
end
if System.verbose == 1
    t = toc;
    disp(['Completed in ' int2str(t) ' seconds !']);
end;
end

