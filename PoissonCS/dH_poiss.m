function [H] = dH_poiss(x)

fun=@(z) (z.*exp(-x.*z)-z)./(z.*log(1-z));
H=-log(x)+integral(@(z)(x>.05).*fun(z),0,1,'ArrayValue',true);
H(isnan(H))=0;
end

