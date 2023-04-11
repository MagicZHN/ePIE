clc;clear;
%重建像素尺寸4.45nm，探针450pixels，步长67pixels，扫描范围603pixels
%初始探针猜测：振幅圆孔，相位平面
m=1630; %样品函数尺寸
n=1630;
X=1:1:m; %大坐标轴
Y=1:1:n;
a=1024; %探针函数盒子尺寸
b=1024;
x=1:1:a; %小坐标轴
y=1:1:b;
aa=0.1; %更新函数参数
bb=0.1;
k=50; %迭代次数
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

probe_gess_all=ones(a,b); %猜测探针函数，全范围
probe_gess_amplitude=probe_gess_all.*circ(a/2,b/2); %猜测探针函数振幅，限制范围
probe_gess_phase=zeros(a,b); %猜测探针函数相位
probe_gess=probe_gess_amplitude.*exp(1i.*probe_gess_phase.*1); %猜测探针函数
% figure(8),imshow(abs(probe_gess));title('猜测探针振幅');%caxis([0,1]);colorbar;colormap(gray);
% figure(9),imshow(angle(probe_gess));title('猜测探针相位');%caxis([-pi,pi]);colorbar;colormap(gray);
%%  重建
for count=1:1:k
    for j=1:1:10
        for i=1:1:10
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

%% 调用函数
function ZZ=circ(i,j)%掩膜生成函数，ij为圆心中心点坐标
m=1024; %矩阵的函数
n=1024; %矩阵的列数
r=225;   %生成圆的半径
ZZ=ones(m,n);
X=1:1:m;
Y=1:1:n;
[XX,YY]=meshgrid(X,Y);
circle=(XX-i).^2+(YY-j).^2; %计算出每一点到圆心的距离的平方
ZZ(circle>r*r)=0;   %找到圆外的元素，并复制为0
ZZ(circle<=r*r)=1;   %找到圆内的元素，并复制为1
end