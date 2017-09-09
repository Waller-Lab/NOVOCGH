function [ amplitude ] = function_fresnelPropField( phase, source, Hstack)
%function_fresnelPropField computes the amplitude determined by the
%Fresnel propagation kernel Hstack.
    objectField = source.*exp(1i * phase);
    amplitude = fftshift(fft2(objectField .* Hstack));
end
 
