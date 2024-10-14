function [obj_output]=calc_obj(S,sup,g_init,itnum,diffraction_type,d,delta,lambda)
if(nargin<=6)
    diffraction_type='fourier';
    d=1;
    delta=0;
    lambda=0;
end

    g=g_init;
    for i = 1:itnum
        G  = propagate(2*g.*sup-g,d,diffraction_type,delta,lambda);
        G2 = S.*exp(1i.*angle(G));
        g2 = propagate(G2,-d,diffraction_type,delta,lambda);
    
        g = g + g2 - g.*sup;
    end
    
    obj_output=g;
    
end

