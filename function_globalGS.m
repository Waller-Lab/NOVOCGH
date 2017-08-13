function [ GS ] = function_globalGS(System, HStacks, Masks )
if System.verbose == 1
    disp('Sequential Gerchberg-Saxton hologram computation begins...');
    tic;
end;
ims = Masks;
[NX,NY,NZ] = size(Masks);
if System.useGPU == 1
    imageInten = zeros(NX,NY,NZ, 'gpuArray');
else
    imageInten = zeros(NX,NY,NZ);
end
[ Superposition] = function_Superoposition( System,HStacks,Masks );
if System.verbose == 1
    fprintf('Completed cycles #');
end;
Superposition.phase = mod(Superposition.phase, 2*pi) - pi;
im = System.source.*exp(1i * Superposition.phase);
for n = 1:System.maxiter+1
    if System.verbose == 1
        if mod(n, 10) == 0
            fprintf(int2str(n));fprintf(',')
        end
    end
    index = 1:NZ;
    tempim = 0 * im;
    for i = index
        imagez = fftshift(fft2(im .* HStacks(:,:,i)));
        target = ims(:,:,i) + System.GSoffset;
        imagez = sqrt(target) .* exp(1i * angle(imagez));
        tempim = tempim + ifft2(ifftshift(imagez))./HStacks(:,:,i);
    end
    im = System.source.*exp(1i * angle(tempim));
end

%phase = gather(angle(im));
if System.verbose == 1
    t = toc;
    disp(['Completed in ' int2str(t) ' seconds !']);
end;

GS.hologram = im;
GS.phase = gather(angle(GS.hologram));
end

