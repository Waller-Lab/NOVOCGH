function [ Superposition ] = function_Superposition( System,HStacks,Masks )

%%% This is the implemetnation of CGH Superposition algorithm. 

if System.verbose == 1; disp('Superposition hologram computation begins...'); tic; end;
[NX,NY,LZ] = size(Masks);
hologram = zeros(NX, NY);
for i = 1:LZ
    target = max(Masks(:,:,i),0);
    rng(50);%MODIFIED
    hologram = hologram + ifft2(ifftshift(sqrt(target) .* exp(1i * randn(NX,NY))))./HStacks(:,:,i);
end

Superposition.hologram = hologram;
Superposition.phase = gather(angle(Superposition.hologram));

if System.verbose == 1
t = toc;
disp(['Superposition Hologram - Completed in ' int2str(t) ' seconds !']);
end
end

