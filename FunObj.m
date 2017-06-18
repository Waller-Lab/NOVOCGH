function [loss, df ] = FunObj( phase, source, z, Nx, Ny, thresholdh, thresholdl, maskFun, fresnelKernelFun, useGPU)
% Computes loss and gradient.
% This function is called by Matlab's fmincon library.



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
    imageInten = abs(imagez.^2);

    
    maskh = mask .* (imageInten < thresholdh);
    maskl = (1-mask) .* (imageInten > thresholdl);
    
    diffh = maskh .* (imageInten - thresholdh);
    diffl = maskl .* (imageInten - thresholdl);
    
    temph = imagez.*diffh;
    temph = conj(HStack).*(Nx*Ny*ifft2(ifftshift(temph)));
    templ = imagez.*diffl;
    templ = conj(HStack).*(Nx*Ny*ifft2(ifftshift(templ)));


    loss = loss + sum(sum(diffh.^2 + diffl.^2)); 
    
    
    df = df +  temph + templ;
end
dfphase = source.*(- real(df).*sin(phase) + imag(df) .* cos(phase));
df = gather(real(dfphase(:)));
loss = gather(real(loss));
end

