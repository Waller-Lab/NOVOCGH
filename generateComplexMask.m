function [ mask ] = generateComplexMask( zi, Nx, Ny, ims, imdepths)
% GENERATEMASK helper function that mathces images to binary target
% intensity.

mask = zeros(Nx, Ny);
for i = 1:numel(imdepths)
    if zi == imdepths(i)
        mask = ims(:,:,i);
        break
    end
end

end

