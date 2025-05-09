function [InversionResult,AfterFitting]=Gaborfitting_ALL(Vecter_data,df,Vecter_p,VecterCut_start,VecterCut_end,Vecter_cp,Vecter_alpha,Vecter_beta,k,EffectiveAtoms_Start,EffectiveAtoms_End,flip)

VecterLength=length(Vecter_data(1,:));
AfterFitting=zeros(size(Vecter_data));
for i=1:VecterLength
    if flip(i)==-1
        [~,AfterFitting(:,i)]=Gaborfitting_EachP_Flip(Vecter_data(:,i),df,Vecter_p(i),VecterCut_start(i),VecterCut_end(i),Vecter_cp(i),Vecter_alpha(i),Vecter_beta(i));
    elseif flip(i)==1
        [~,AfterFitting(:,i)]=Gaborfitting_EachP(Vecter_data(:,i),df,Vecter_p(i),VecterCut_start(i),VecterCut_end(i),Vecter_cp(i),Vecter_alpha(i),Vecter_beta(i));
    end  
end
G=AfterFitting;

InversionResult=zeros(size(G(:,1)));
for i=1:EffectiveAtoms_End-EffectiveAtoms_Start+1
    InversionResult=G(:,i).*k(i)+InversionResult;
end