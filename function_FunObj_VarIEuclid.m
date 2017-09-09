function [loss, df ] = function_FunObj_VarIEuclid( phase, source, z, Nx, Ny,  maskFun, fresnelKernelFun, useGPU)
% 
% This function is called by Matlab's fmincon library.
% Computes loss and gradient for NOVOCGH with a Euclidean cost function, accepts variable intensity patterns.. 

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
    
    I = abs(imagez.^2);
    mass = sum(I(:));
    I = I/mass;
    if sum(mask(:))>0
    V = mask/sum(mask(:)); 
    else 
        V = mask;
    end
    diffh = (I-V).^2;
    temph = Nx*Ny*mass*(I-V);   
    temph = imagez.*temph;
    temph = conj(HStack).*(Nx*Ny*ifft2(ifftshift(temph)));
        %Compute losses
    loss = loss + sum(sum(diffh )); 
        %Compute gradient
    df = df +  temph ;
end
dfphase = source.*(- real(df).*sin(phase) + imag(df) .* cos(phase));
df = gather(real(dfphase(:)));
loss = gather(real(loss));
end

