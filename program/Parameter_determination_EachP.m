function [cp,alpha,beta]=Parameter_determination_EachP(df,amplitude_spectrum_seismic,p,mu)
n=length(amplitude_spectrum_seismic);
S_new=abs(amplitude_spectrum_seismic).^p+0.0;%保证log有意义
S0=S_new./sum(S_new)*df;
F=zeros(n,1);
F(1)=S0(1);
for i=2:n
    F(i)=F(i-1)+S0(i)*df;
end
S0=abs(amplitude_spectrum_seismic);
S=log(S0);
A1=ones(n,1);A2=log(F);A3=log(1-F);
A=[A1,A2,A3];
m=inv(A'*A+mu*eye(3))*A'*S;
cp=exp(m(1));alpha=m(2);beta=m(3);
