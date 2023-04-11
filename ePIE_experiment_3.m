clc;clear;
%�ؽ����سߴ�4.45nm��̽��450pixels������67pixels��ɨ�跶Χ603pixels
%��ʼ̽��²⣺�������Ƭģ�⣬��λ����
m=1630; %��Ʒ�����ߴ�
n=1630;
X=1:1:m; %��������
Y=1:1:n;
a=1024; %̽�뺯�����ӳߴ�
b=1024;
x=1:1:a; %С������
y=1:1:b;
aa=0.1; %���º�������
bb=0.01;
k=50; %��������
%% ̽������
lamda=1.77e-9;                            %���ò���Ϊ����m
D=200e-6;                                  %���÷���������Ƭֱ��Ϊ����m
f=3.39e-3;                                   %���ý���()m
dx=35.6e-9;                                 %һ�����ص��С()m
N=5620;                                   %�����С
kk=2*pi/lamda;                             %��ʸ
dz=33.87e-6;                                %�뽹����()m
z=f+dz;                                   %��������
aaa=128;
% xxx=1:1:aaa;
[xx,yy]=meshgrid(linspace(-N*dx/2,N*dx/2,N)); %�ռ�����
%% ��λ����
dxx=4.45e-9;
[xx2,yy2]=meshgrid(linspace(-a*dxx/2,a*dxx/2,a)); %�ռ�����
A1=1./sqrt(xx2.^2+yy2.^2+dz.^2).*exp(1j*kk.*sqrt(xx2.^2+yy2.^2+dz.^2));
phase=angle(A1);
%% ���Ʋ���Ƭ
M=zeros(N,N);
for mm=1:N    
    for nn=1:N             %���ö���ѭ�����������������Ƭ��Ļ�ϸ���
        p=sqrt(xx(1,mm).^2+yy(nn,1).^2);        %���������Բ�뾶
        kn=fix(p.^2./(lamda.*f));             %�����벨����
        if p>D/2                            %�����Ļ�ϵ���ڲ���Ƭ�뾶
            M(mm,nn)=0;                    % ���ú�ɫ��ʾ���������ǲ���Ƭ����
        else
            if mod(kn,2)==1                %�жϰ벨������ż��Ϊ������Ϳ��
                M(mm,nn)=0;
            else
                M(mm,nn)=1;                      %Ϊż����Ϳ�ڱ�ʾ͸��
            end
        end
    end
end
% figure(1),imshow(M);                                  %���Ʋ���Ƭͼ��
%% ���ιⳡ
h=exp(1j*kk*z).*exp((1j*kk*(xx.^2+yy.^2))./(2*z))./(1j*lamda*z);
q=exp((1j*kk*(xx.^2+yy.^2))./(2*z));
A2=h.*fftshift(fft2(M.*q));
amplitude_all=abs(A2)./max(max(abs(A2)));
%% ��ʾ���ľֲ�
amplitude=amplitude_all((N-aaa)/2+1:aaa+(N-aaa)/2,(N-aaa)/2+1:aaa+(N-aaa)/2);
%% ��չ
amplitude_nice=imresize(amplitude,[1024,1024],'bicubic');
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

probe_gess_amplitude=amplitude_nice; %�²�̽�뺯�����
probe_gess_phase=phase; %�²�̽�뺯����λ
probe_gess=probe_gess_amplitude.*exp(1i.*probe_gess_phase.*1); %�²�̽�뺯��
% figure(8),imshow(abs(probe_gess));title('�²�̽�����');%caxis([0,1]);colorbar;colormap(gray);
% figure(9),imshow(angle(probe_gess));title('�²�̽����λ');%caxis([-pi,pi]);colorbar;colormap(gray);
%%  �ؽ�
for count=1:1:k
    for j=10:-1:1
        for i=10:-1:1
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