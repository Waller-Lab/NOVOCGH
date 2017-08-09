function [ amplitude ] = function_fresnelPropField( phase, source, Hstack)
%FRESNELPROP Compute the generated amplitude determined by the
%Fresnel propagation kernel Hstack.
    objectField = source.*exp(1i * phase);
    amplitude = fftshift(fft2(objectField .* Hstack));
end
 
