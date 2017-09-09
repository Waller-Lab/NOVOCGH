function [loss, df ] = function_FunObj_VarI( phase, source, z, Nx, Ny, thresholdh, thresholdl, maskFun, fresnelKernelFun, useGPU)
% 
% This function is called by Matlab's fmincon library.
% Computes loss and gradient for NOVOCGH with a threshold-based cost function, with variable intensity patterns, and enforces darkness wherever specified ! 

if useGPU
df = zeros(Nx, Ny, 'gpuArray');
phase = gpuArray(phase);
source = gpuArray(source);
else
df = zeros(Nx, Ny);
end

loss = 0; 
phase = reshape(phase, [Nx, Ny]);
objectField = source.*exp(1i * phase);

for i = 1 : numel(z)
    HStack = fresnelKernelFun(i,i);
    mask = maskFun(i,i);  
   
    imagez = fftshift(fft2(objectField .* HStack));
    intenmask = mask;
    intenmask(intenmask==0) = 1; 
    %imageInten = abs(imagez.^2)./intenmask;
   imageInten = (abs(imagez.^2)./intenmask).*(mask >0); 
    maskh = (mask>0) .* (imageInten < thresholdh);
    diffh = maskh .* (imageInten - thresholdh);
    temph = imagez.*diffh;
    temph = conj(HStack).*(Nx*Ny*ifft2(ifftshift(temph)));
        %Compute losses
    loss = loss + sum(sum(diffh.^2 )); 
        %Compute gradient
    df = df +  temph;
end
dfphase = source.*(- real(df).*sin(phase) + imag(df) .* cos(phase));
df = gather(real(dfphase(:)));
loss = gather(real(loss));
end

