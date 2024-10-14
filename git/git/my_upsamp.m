function new_obj=my_upsamp(obj,sup_d,N_d,m_d,N,m)
new_obj=obj.*sup_d;
new_obj=new_obj((N_d-m_d)/2+1:(N_d+m_d)/2,(N_d-m_d)/2+1:(N_d+m_d)/2);
new_obj=padarray(imresize(new_obj,[m,m]),[(N-m)/2,(N-m)/2]);
end