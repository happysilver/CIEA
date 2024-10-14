close all;clear;clc
addpath(genpath('functions'))

diffraction_type='fourier';
d=1;
opti_method="ga";

N = 1024;
m = 256;

sample_p = double(imread('data\images\Cameraman.tif'));
sample_p = sample_p(:,:,1);
sample_p = padarray(imresize(sample_p,[m,m]),[(N-m)/2,(N-m)/2],0);
sample_p = MatMap(sample_p,0,1);
sample=exp(1i.*sample_p);
sup=triMask(N,1.5*m,N/2,N/2);
I0=abs(propagate(sample.*sup,1));
I0=MatMap(I0,0,1);

ratio=8;
block_size=8;
samp_num=3;
N_d=N/ratio;
m_d=floor(m/ratio);

T=dctmtx(block_size);
sup_d = imresize(sup,[N_d,N_d]);

x_num=(m_d/block_size)^2*samp_num*2;

block_num=m_d/block_size;
index=0;
for ii=1:block_num
    for jj=1:block_num
        index=index+1;x_row(index)=(ii-1)*block_size+1;x_col(index)=(jj-1)*block_size+1;
        index=index+1;x_row(index)=(ii-1)*block_size+1;x_col(index)=(jj-1)*block_size+2;
        index=index+1;x_row(index)=(ii-1)*block_size+2;x_col(index)=(jj-1)*block_size+1;
    end
end
x_index=sub2ind([m_d,m_d],x_row,x_col);

x_0=ones(1,x_num);
x_lb=ones(size(x_0))*(-3);
x_ub=ones(size(x_0))*(7);

I1=I0((N-N_d)/2+1:(N+N_d)/2,(N-N_d)/2+1:(N+N_d)/2);
I1=I1./mean2(I1)/2;
least_func=@(x) least_func_whole(x,sup_d,d,diffraction_type,I1,block_size,N_d,m_d,x_index,T);
x_out=x_0;

if(opti_method=="fmincon")
    options = optimoptions('fmincon','MaxFunctionEvaluations',1e4,'MaxIterations',1e5,'ConstraintTolerance',1e-8,'StepTolerance',1e-12);
    x_out=fmincon(@least_func,x_0,[],[],[],[],x_lb,x_ub,[],options);
elseif(opti_method=="ga")
    options = optimoptions('ga','PlotFcn',@gaplotbestf,'FunctionTolerance',1e-5);
    x_out=ga(least_func,x_num,[],[],[],[],x_lb,x_ub,[],options);
elseif(opti_method=="PSO")
    options = optimoptions('particleswarm','PlotFcn','pswplotbestf');
    x_out=particleswarm(least_func,x_num,x_lb,x_ub,options);
elseif(opti_method=="simulannealbnd")
    options = optimoptions('simulannealbnd','PlotFcn','saplotbestf');
    x_out=simulannealbnd(least_func,x_0,x_lb,x_ub,options);    
elseif(opti_method=="EGO")
    %
end

oo=make_obj_whole(x_out,block_size,N_d,m_d,x_index,T).*sup_d;
disp('down opti ok')

I1=I0((N-N_d)/2+1:(N+N_d)/2,(N-N_d)/2+1:(N+N_d)/2);
I1=MatMap(I1,0,1);
itnum=50;
[obj_out0]=calc_obj(I1,sup_d,make_obj_whole(x_out,block_size,N_d,m_d,x_index,T).*sup_d,itnum,diffraction_type);
new_obj_0=my_upsamp(obj_out0,sup_d,N_d,m_d,N,m);
obj_show=(N-m)/2+1:(N+m)/2;
figure();imshow(angle(new_obj_0(obj_show,obj_show)).*sup(obj_show,obj_show),[]),title('CIEA init')
disp('additional init process ok')

I1=I0;
itnum=50;
[obj_out1]=calc_obj(I1,sup,exp(1i.*rand(N,N)),itnum,diffraction_type);
[obj_out2]=calc_obj(I1,sup,new_obj_0,itnum,diffraction_type);
obj_show=(N-m)/2+1:(N+m)/2;
figure();
subplot(1,2,1),imshow(angle(obj_out1(obj_show,obj_show)).*sup(obj_show,obj_show),[]),title('AP')
subplot(1,2,2),imshow(angle(obj_out2(obj_show,obj_show)).*sup(obj_show,obj_show),[]),title('CIEA')
disp('end')

function y=least_func_whole(x,sup_d,d,diffraction_type,I1,block_size,N_d,m_d,x_index,T)
I_temp=abs(propagate(make_obj_whole(x,block_size,N_d,m_d,x_index,T).*sup_d,d,diffraction_type));
I_temp=I_temp./mean2(I_temp)/2;
y=norm(I_temp-I1,2);
end