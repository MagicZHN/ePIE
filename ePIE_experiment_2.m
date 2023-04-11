clc;clear;
%�ؽ����سߴ�4.45nm��̽��450pixels������67pixels��ɨ�跶Χ603pixels
%��ʼ̽��²⣺���Բ�ף���λƽ��
m=1630; %��Ʒ�����ߴ�
n=1630;
X=1:1:m; %��������
Y=1:1:n;
a=1024; %̽�뺯�����ӳߴ�
b=1024;
x=1:1:a; %С������
y=1:1:b;
aa=0.1; %���º�������
bb=0.1;
k=50; %��������
%% ��ȡ����
I=zeros(a,b,10*10);
load 700eV2um300nm220ms10x10_5.mat
for ii=1:1:100
    I(:,:,ii)=double(ccdfs(:,:,ii)); %��ȡ����ͼ
end
%% ��ʼ���²�
object_gess_amplitude=ones(m,n).*0.36; %��ʼ�Ĳ²���Ʒ�������
object_gess_phase=ones(m,n).*(-0.4); %��ʼ�Ĳ²���Ʒ������λ
object_gess=object_gess_amplitude.*exp(1i.*object_gess_phase); %��ʼ�Ĳ²���Ʒ����

probe_gess_all=ones(a,b); %�²�̽�뺯����ȫ��Χ
probe_gess_amplitude=probe_gess_all.*circ(a/2,b/2); %�²�̽�뺯����������Ʒ�Χ
probe_gess_phase=zeros(a,b); %�²�̽�뺯����λ
probe_gess=probe_gess_amplitude.*exp(1i.*probe_gess_phase.*1); %�²�̽�뺯��
% figure(8),imshow(abs(probe_gess));title('�²�̽�����');%caxis([0,1]);colorbar;colormap(gray);
% figure(9),imshow(angle(probe_gess));title('�²�̽����λ');%caxis([-pi,pi]);colorbar;colormap(gray);
%%  �ؽ�
for count=1:1:k
    for j=1:1:10
        for i=1:1:10
            s=10*(j-1)+i;  %��ǰɨ��λ��������
            object_gess_r=object_gess(67*i-66:a+67*i-67,67*j-66:a+67*j-67); %�Ӳ²���Ʒ������ȡ����һ��
            object_gess_r_save=object_gess_r; %�����ֲ��²���Ʒ
            probe_gess_save=probe_gess; %�����²�̽�뺯��
            o_max=max(max((abs(object_gess_r)).^2));
            p_max=max(max((abs(probe_gess)).^2));
            
            phy_gess=object_gess_r.*probe_gess; %�²���䲨
            phyF_gess=fftshift(fft2(phy_gess)); %�²����䲨
            phyF_update=sqrt(I(:,:,s)).*phyF_gess./abs(phyF_gess); %����滻��������λ
            phy_update=ifft2(fftshift(phyF_update)); %�渵��Ҷ�任��ʵ�ռ�
            
            object_gess_r=object_gess_r+(conj(probe_gess_save)./p_max).*aa.*(phy_update-phy_gess); %����
            probe_gess=probe_gess+(conj(object_gess_r_save)./o_max).*bb.*(phy_update-phy_gess);

            object_gess(67*i-66:a+67*i-67,67*j-66:a+67*j-67)=object_gess_r;%�ֲ���Ʒ������ԭ�ز²���Ʒ����
        end
    end
end
%% ��ʾ
figure(1),imshow(abs(object_gess));title('�ؽ���Ʒ���');%caxis([0,1]);colorbar;colormap(gray);
% figure(6),mesh(X,Y,abs(object_gess));title('�ؽ���Ʒ���')
figure(2),imshow(angle(object_gess));title('�ؽ���Ʒ��λ');%caxis([-pi,pi]);colorbar;colormap(gray);
% figure(3),mesh(x,y,abs(probe_gess));title('�ؽ�̽�����')
figure(4),imshow(abs(probe_gess));title('�ؽ�̽�����');%caxis([0,1]);colorbar;colormap(gray);
% figure(5),mesh(x,y,angle(probe_gess));title('�ؽ�̽����λ')
figure(7),imshow(angle(probe_gess));title('�ؽ�̽����λ')

%% ���ú���
function ZZ=circ(i,j)%��Ĥ���ɺ�����ijΪԲ�����ĵ�����
m=1024; %����ĺ���
n=1024; %���������
r=225;   %����Բ�İ뾶
ZZ=ones(m,n);
X=1:1:m;
Y=1:1:n;
[XX,YY]=meshgrid(X,Y);
circle=(XX-i).^2+(YY-j).^2; %�����ÿһ�㵽Բ�ĵľ����ƽ��
ZZ(circle>r*r)=0;   %�ҵ�Բ���Ԫ�أ�������Ϊ0
ZZ(circle<=r*r)=1;   %�ҵ�Բ�ڵ�Ԫ�أ�������Ϊ1
end