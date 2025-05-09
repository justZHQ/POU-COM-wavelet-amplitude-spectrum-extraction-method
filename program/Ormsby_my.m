function [t,ormsby,f,amplitude_spectrum]=Ormsby_my(dt,L,F_L1,F_L2,F_H1,F_H2,fa,fmax)
% dt;%ʱ�����/s
% L;%�Ӳ����ȣ�%wavelength=2*L+1
% F_L1,F_L2,F_H1,F_H2;%ormsby�Ӳ���Ƶ��/HZ������F_L1<F_L2<F_H1<F_H2
% fa;%��λ
% fmax;%�����ʾƵ��/Hz
% t;%���ʱ������/s
% ormsby;%���ʱ����ormsby�Ӳ�
% f;%���Ƶ������
% amplitude_spectrum;%���ormsby�Ӳ��������

t=-L*dt:dt:L*dt;
t((length(t)+1)/2)=0.000000001;
s1=(((pi*F_H2).^2/(pi*F_H2-pi*F_H1)).*((sin(pi*F_H2.*t)./(pi*F_H2.*t)).^2)-((pi*F_H1).^2/(pi*F_H2-pi*F_H1)).*((sin(pi*F_H1.*t)./(pi*F_H1.*t)).^2))...
    -(((pi*F_L2).^2/(pi*F_L2-pi*F_L1)).*((sin(pi*F_L2.*t)./(pi*F_L2.*t)).^2)-((pi*F_L1).^2/(pi*F_L2-pi*F_L1)).*((sin(pi*F_L1.*t)./(pi*F_L1.*t)).^2));
s1=s1/max(s1);


s2=hilbert(s1);
ormsby=real(s2)*cos(fa)+imag(s2)*sin(fa);
[f,amplitude_spectrum]=Amplitude_spectrum_my(dt,ormsby,fmax);

% figure;
% plot(real(s1));
% figure;
% plot(amplitude_spectrum);
% 
% figure;