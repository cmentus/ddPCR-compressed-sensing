function [x_estimate] = decode_sparse(Pool_values,Pool_mat,alpha,num_iter,support,N)
NN=length(support);
tau=zeros(NN,1);
Pool_mat=Pool_mat(:,support)+.0001;
cvx_begin
    expression x_e(NN)
    variable x_n(NN)
    minimize(sum(Pool_mat*x_n-Pool_values.*log(Pool_mat*x_n))+alpha*sum(tau.*x_n));
    subject to
        x_n>=ones(NN,1)*.0001
cvx_end
x_e=x_n;
for iter=1:num_iter-1

    tau=.5*x_e.^(-.5);
    cvx_begin
        variable x_n(NN)
        minimize( sum(Pool_mat*x_n-Pool_values.*log(Pool_mat*x_n))...
            + alpha*sum(tau.*x_n) );
        subject to
            x_n >= 0.0001*ones(NN,1);
    cvx_end
    x_e=x_n;
end
x_estimate=zeros(N,1);
x_estimate(support)=x_n;

end

