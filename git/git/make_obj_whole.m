function obj_out=make_obj_whole(x,block_size,N_d,m_d,x_index,T)
    obj_samp_angle=zeros(m_d,m_d);
    obj_samp_angle(x_index)=x(1:length(x)/2);
    obj_samp_abs=zeros(m_d,m_d);
    obj_samp_abs(x_index)=x(length(x)/2+1:end);
    fun_DCT=@(x) T'*x.data*T;
    obj_out_angle=blockproc(obj_samp_angle,[block_size,block_size],fun_DCT);
    obj_out_abs=blockproc(obj_samp_abs,[block_size,block_size],fun_DCT);
    obj_out=obj_out_abs.*exp(1i.*obj_out_angle);

    obj_out=padarray(obj_out,[(N_d-m_d)/2,(N_d-m_d)/2],0);
end