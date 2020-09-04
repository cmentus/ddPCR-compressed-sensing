function [x_n] = decode_l1(Pool_values,Pool_mat,alpha,support,N)
% use support
Pool_mat=Pool_mat(:,support)+.00001;
NN=length(support);
cvx_begin
    expression x_e(N)
    variable x_estimate(NN)
    minimize(sum(Pool_mat*x_estimate-Pool_values.*log(Pool_mat*x_estimate))...
        + alpha * sum(x_estimate))
    subject to
        x_estimate>=0.0001;
cvx_end
x_n=zeros(N,1);
x_n(support)=x_estimate;

end

