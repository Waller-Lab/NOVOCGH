function [ Intensity ] = function_VolumeIntensity( System,phase,HStacks)
if System.verbose == 1
tic;
fprintf('Computing volume image...'); end;
[NX,NY,NZ] = size(HStacks);
Intensity = zeros(NX,NY,NZ);
phase = mod(phase, 2*pi) - pi;
for i = 1:NZ
    imagez = function_fresnelProp(phase, System.source, HStacks(:,:,i));
    Intensity(:,:,i) = gather(imagez);
end
%Intensity = Intensity/max(Intensity(:));
if System.verbose == 1
t = toc;
disp(['Completed in ' int2str(t) ' seconds !']); 
end;
end

