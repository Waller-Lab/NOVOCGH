function [ Scores ]= function_ScoresI(KickMasks,Masks, IntensityStack);
A = IntensityStack.*(Masks>0);
lambda = sum(sum(sum(A.*Masks)))/sum(sum(sum(A.*A)));
A = lambda*A;

Scores.NegativeError = squeeze(gather(sum(sum(IntensityStack.*(KickMasks>0),1),2)./sum(sum(IntensityStack,1),2))); %from 0 to 1, amount of misplaced light
Scores.PositiveError = squeeze(gather(sum(sum(IntensityStack.*(Masks>0),1),2)./sum(sum(IntensityStack,1),2))); %from 0 to 1, amount of light placed in targets light
Scores.TargetMatcherror = squeeze(gather(sum(sum(abs(A-Masks),1),2)./sum(sum(A,1),2)));

Scores.NegativeError(isnan(Scores.NegativeError)) = 0;
Scores.PositiveError(isnan(Scores.PositiveError)) = 0;
Scores.TargetMatcherror(isnan(Scores.TargetMatcherror)) = 0;


end

