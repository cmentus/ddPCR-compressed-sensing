function [H] = H_poiss(x)

fun=@(z) (1-exp(-x.*z)-x.*z)./(z.*log(1-z));
H=-x.*log(x./exp(1))+integral(@(z)(x>.05).*fun(z),0,1,'ArrayValue',true);
H(isnan(H))=0;
end

