function [GaborIntransformation]=InGabor(EffectiveAtoms_Start,EffectiveAtoms_End,GaborTimeShiftedSignalTruncation)

%%
G=zeros(length(GaborTimeShiftedSignalTruncation(:,1)),EffectiveAtoms_End-EffectiveAtoms_Start+1);
GaborIntransformation=zeros(size(G(:,1)));
for j=EffectiveAtoms_Start:EffectiveAtoms_End
    GaborIntransformation=GaborIntransformation+GaborTimeShiftedSignalTruncation(:,j);    
end
