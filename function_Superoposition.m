function [ hologram ] = function_Superoposition( System,HStacks,Masks )
%This algorithm conducts a simple superposition 
if System.verbose == 1; fprintf('Superposition hologram computation begins...'); tic; end;
[NX,NY,LZ] = size(Masks);
hologram = zeros(NX, NY);
for i = 1:LZ
    target = max(Masks(:,:,i),0);
    hologram = hologram + ifft2(ifftshift(sqrt(target) .* exp(1i * randn(NX,NY))))./HStacks(:,:,i);
end
if System.verbose == 1
t = toc;
disp(['Completed in ' int2str(t) ' seconds !']);
end
end

