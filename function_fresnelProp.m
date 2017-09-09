function [ imagez ] = function_fresnelProp( phase, source, Hstack)
%function_fresnelProp computes the generated intensity image determined by the Fresnel propagation kernel Hstack.
    objectField = source.*exp(1i * phase);
    imagez = abs(fftshift(fft2(objectField .* Hstack)).^2);
end
 
