function I_out=propagate(I_in,d,type,delta,lambda)
if(nargin<3)
    type='fourier';
end

switch type
    case 'fourier'
        if(d>0)
            I_out=fftshift(fft2(I_in));
        else
            I_out=ifft2(ifftshift(I_in));
        end
    case 'fresnel'
        I_out=fresnel_advance(I_in, delta, delta, d, lambda);
    otherwise
 
end
end