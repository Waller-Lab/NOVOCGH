function [HStack]=GenerateFresnelPropagationStack(Nx,Ny,z, lambda, psx, psy, useGPU)
% Lambda, ps, z has unit meter.
% z is the distance from focal plane.
% lambda is wavelength
% ps is pixel size.

cx=[1:Nx] - (floor(Nx/2)+1);
cy=[1:Ny] - (floor(Ny/2)+1);


if useGPU
    cx = gpuArray(cx); cy = gpuArray(cy);
end

[us, vs]=ndgrid(cx, cy);
us=us/Nx/psx; vs=vs/Ny/psy;

HStack = exp(1i*lambda*pi*(us.^2+vs.^2)* z/2);

end

