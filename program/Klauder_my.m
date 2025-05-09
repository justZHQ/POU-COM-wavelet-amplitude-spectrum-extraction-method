function [t,klauder,f,amplitude_spectrum]=Klauder_my(dt,L,F_L,F_H,T,fa,fmax)
% dt;%ʱ�����/s
% L;%�Ӳ����ȣ�%wavelength=2*L+1
% F_L,F_H;%klauder�Ӳ��ĵͽ�Ƶ�͸߽�Ƶ/HZ������F_L<F_H
% T;%����/s
% fa;%��λ
% fmax;%�����ʾƵ��/Hz
% t;%���ʱ������/s
% klauder;%���ʱ����klauder�Ӳ�
% f;%���Ƶ������
% amplitude_spectrum;%���klauder�Ӳ��������
t=-L*dt:dt:L*dt;
t((length(t)+1)/2)=0.000000001;
k=(F_H-F_L)/T;
f0=(F_H+F_L)/2;
s1=sin(pi*k.*t.*(T-t))./(pi*k.*t).*cos(2*pi*f0.*t);
s1=s1/max(s1);
s2=hilbert(s1);
klauder=real(s2)*cos(fa)+imag(s2)*sin(fa);
[f,amplitude_spectrum]=Amplitude_spectrum_my(dt,klauder,fmax);
