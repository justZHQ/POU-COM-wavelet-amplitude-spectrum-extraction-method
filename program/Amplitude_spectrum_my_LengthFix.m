function [f,F,amplitude_spectrum,N]=Amplitude_spectrum_my_LengthFix(dt,signal,fmax,number)
% dt;%ʱ�����/s
% signal;%����ʱ���ź�
% fmax;%�����ʾƵ��/Hz
% f;%���Ƶ������/Hz
% amplitude_spectrum;%�����ź�signal�������
F=fft(signal,2^ceil(log2(number)));
N=length(F);
amplitude_spectrum=abs(F);%abs(Cn)=abs(X(k))=sqrt(realX(k)^2+imagX(k)^2)
f=(0:N-1)*(1/(N*dt));%Ƶ�ʲ���������ڻ���Ƶ��f0=1/T;T=dt*N=N/fs
fmax_number=ceil(fmax/(1/(N*dt)));
f=f(1:fmax_number);
amplitude_spectrum=amplitude_spectrum(1:fmax_number)*2/N;%Ƶ������ʵ���֮��Ĺ�ϵ��������ף�
F=F(1:fmax_number);
end
