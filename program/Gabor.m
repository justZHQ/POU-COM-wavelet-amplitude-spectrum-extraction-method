function [NormalizedGaborAtomicLibrary,TimeShiftedSignal,GaborTimeShiftedSignalTruncation]=Gabor(dt,GaussianLength,GausiannVariance,GivenSignalLength,Scale,TimeOffset,Signal)
T=-GaussianLength*dt:dt:dt*GaussianLength;
Gaussian=exp(-T.^2/(2*GausiannVariance^2));

GaborAtomicLibrary=zeros(length(Gaussian)+GivenSignalLength-1,floor(GivenSignalLength/Scale));
k=1;
for i=1:Scale:GivenSignalLength
    GaborAtomicLibrary(i:length(Gaussian)+i-1,k)=Gaussian;
    k=k+1;
end
AtomicLibraryNumber=k-1;
%% 
SumAtomic=sum(GaborAtomicLibrary,2);%按列求和
NormalizedGaborAtomicLibrary=zeros(size(GaborAtomicLibrary));
for i=1:AtomicLibraryNumber
    NormalizedGaborAtomicLibrary(:,i)=GaborAtomicLibrary(:,i)./SumAtomic;
end
%%
OutputSignalLength=length(GaborAtomicLibrary(:,1));
TimeShiftedSignal=zeros(OutputSignalLength,1);
TimeShiftedSignal(TimeOffset:length(Signal)+TimeOffset-1)=Signal;
GaborTimeShiftedSignal=zeros(length(Signal),AtomicLibraryNumber);
for i=1:AtomicLibraryNumber
    GaborTimeShiftedSignal(:,i)=NormalizedGaborAtomicLibrary(TimeOffset+1:TimeOffset+length(Signal),i).*Signal;
end
%%
GaborTimeShiftedSignalTruncation=GaborTimeShiftedSignal;
end
