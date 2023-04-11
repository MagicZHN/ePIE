clc;clear;
%重建像素尺寸4.45nm，探针450pixels，步长67pixels，扫描范围603pixels
%初始探针猜测：振幅波带片模拟，相位球面
m=1630; %样品函数尺寸
n=1630;
X=1:1:m; %大坐标轴
Y=1:1:n;
a=1024; %探针函数盒子尺寸
b=1024;
x=1:1:a; %小坐标轴
y=1:1:b;
aa=0.1; %更新函数参数
bb=0.01;
k=50; %迭代次数
%% 探针生成
lamda=1.77e-9;                            %设置波长为（）m
D=200e-6;                                  %设置菲涅尔波带片直径为（）m
f=3.39e-3;                                   %设置焦距()m
dx=35.6e-9;                                 %一个像素点大小()m
N=5620;                                   %矩阵大小
kk=2*pi/lamda;                             %波矢
dz=33.87e-6;                                %离焦距离()m
z=f+dz;                                   %传播距离
aaa=128;
% xxx=1:1:aaa;
[xx,yy]=meshgrid(linspace(-N*dx/2,N*dx/2,N)); %空间坐标
%% 相位生成
dxx=4.45e-9;
[xx2,yy2]=meshgrid(linspace(-a*dxx/2,a*dxx/2,a)); %空间坐标
A1=1./sqrt(xx2.^2+yy2.^2+dz.^2).*exp(1j*kk.*sqrt(xx2.^2+yy2.^2+dz.^2));
phase=angle(A1);
%% 绘制波带片
M=zeros(N,N);
for mm=1:N    
    for nn=1:N             %设置二重循环依次求菲涅尔波带片屏幕上各点
        p=sqrt(xx(1,mm).^2+yy(nn,1).^2);        %求各点所在圆半径
        kn=fix(p.^2./(lamda.*f));             %求各点半波点数
        if p>D/2                            %如果屏幕上点大于波带片半径
            M(mm,nn)=0;                    % 则用黑色表示背景，不是波带片部分
        else
            if mod(kn,2)==1                %判断半波带数奇偶，为奇数则涂黑
                M(mm,nn)=0;
            else
                M(mm,nn)=1;                      %为偶数则不涂黑表示透光
            end
        end
    end
end
% figure(1),imshow(M);                                  %绘制波带片图像
%% 下游光场
h=exp(1j*kk*z).*exp((1j*kk*(xx.^2+yy.^2))./(2*z))./(1j*lamda*z);
q=exp((1j*kk*(xx.^2+yy.^2))./(2*z));
A2=h.*fftshift(fft2(M.*q));
amplitude_all=abs(A2)./max(max(abs(A2)));
%% 显示中心局部
amplitude=amplitude_all((N-aaa)/2+1:aaa+(N-aaa)/2,(N-aaa)/2+1:aaa+(N-aaa)/2);
%% 扩展
amplitude_nice=imresize(amplitude,[1024,1024],'bicubic');
%% 读取数据
I=zeros(a,b,10*10);
load 700eV2um300nm220ms10x10_5.mat
for ii=1:1:100
    I(:,:,ii)=double(ccdfs(:,:,ii)); %读取衍射图
end
%% 初始化猜测
object_gess_amplitude=ones(m,n).*0.36; %初始的猜测样品函数振幅
object_gess_phase=ones(m,n).*(-0.4); %初始的猜测样品函数相位
object_gess=object_gess_amplitude.*exp(1i.*object_gess_phase); %初始的猜测样品函数

probe_gess_amplitude=amplitude_nice; %猜测探针函数振幅
probe_gess_phase=phase; %猜测探针函数相位
probe_gess=probe_gess_amplitude.*exp(1i.*probe_gess_phase.*1); %猜测探针函数
% figure(8),imshow(abs(probe_gess));title('猜测探针振幅');%caxis([0,1]);colorbar;colormap(gray);
% figure(9),imshow(angle(probe_gess));title('猜测探针相位');%caxis([-pi,pi]);colorbar;colormap(gray);
%%  重建
for count=1:1:k
    for j=10:-1:1
        for i=10:-1:1
            s=10*(j-1)+i;  %当前扫描位置索引号
            object_gess_r=object_gess(67*i-66:a+67*i-67,67*j-66:a+67*j-67); %从猜测样品函数中取出第一块
            object_gess_r_save=object_gess_r; %保留局部猜测样品
            probe_gess_save=probe_gess; %保留猜测探针函数
            o_max=max(max((abs(object_gess_r)).^2));
            p_max=max(max((abs(probe_gess)).^2));
            
            phy_gess=object_gess_r.*probe_gess; %猜测出射波
            phyF_gess=fftshift(fft2(phy_gess)); %猜测衍射波
            phyF_update=sqrt(I(:,:,s)).*phyF_gess./abs(phyF_gess); %振幅替换，保留相位
            phy_update=ifft2(fftshift(phyF_update)); %逆傅里叶变换回实空间
            
            object_gess_r=object_gess_r+(conj(probe_gess_save)./p_max).*aa.*(phy_update-phy_gess); %更新
            probe_gess=probe_gess+(conj(object_gess_r_save)./o_max).*bb.*(phy_update-phy_gess);

            object_gess(67*i-66:a+67*i-67,67*j-66:a+67*j-67)=object_gess_r;%局部样品函数复原回猜测样品函数
        end
    end
end
%% 显示
figure(1),imshow(abs(object_gess));title('重建样品振幅');%caxis([0,1]);colorbar;colormap(gray);
% figure(6),mesh(X,Y,abs(object_gess));title('重建样品振幅')
figure(2),imshow(angle(object_gess));title('重建样品相位');%caxis([-pi,pi]);colorbar;colormap(gray);
% figure(3),mesh(x,y,abs(probe_gess));title('重建探针振幅')
figure(4),imshow(abs(probe_gess));title('重建探针振幅');%caxis([0,1]);colorbar;colormap(gray);
% figure(5),mesh(x,y,angle(probe_gess));title('重建探针相位')
figure(7),imshow(angle(probe_gess));title('重建探针相位')