

N=238; M=48; % N is the population size, % M is the number of groups
Phi_init=(rand(M,N)<.25)*.99; % Initialization of the optimization algorithm
Phi_init=bsxfun(@rdivide, Phi_init,sqrt(sum(Phi_init.^2,1))); % l2 Normalize the columns
%%
% Start with Phi_init and optimize it. Keep track of average and mutual coherences
Phi=Phi_init;
coherences=[mutual_coherence(Phi)];
avg_coherences=[average_coherence(Phi)]
n_iter=20;

for iteration=1:n_iter
    funs=[];
    Phi_old=Phi;
    for n=1:N
        Gram_matrix=transpose(Phi)*Phi;
        T_n=sqrt(1-max(triu(Gram_matrix,1).^2+tril(Gram_matrix,-1).^2,[],2));
        mu_n=@(f) coherence_with_pool(f,Phi,n);
        f_init=Phi(:,n);
        constrnl=@(f) constr_ball(f,Phi(:,n),T_n(n));
        
        cvx_begin
            variable f_n(M)
            minimize(mu_n(f_n))
            subject to
                norm(f_n-f_init)<=T_n(n);
                f_n>=0;
                %sum(Phi-Phi(:,n),2)+f_n=ones(M,1)
        cvx_end
        T_n(n);
        Phi(:,n)=f_n/norm(f_n);
        funs=[funs [mu_n(f_init); mu_n(f_n)]];
        coherences=[coherences mutual_coherence(Phi)];
        avg_coherences=[avg_coherences average_coherence(Phi)];
        
    end
    
    %coherences=[coherences mutual_coherence(Phi_new)];
    %plot(coherences)
end





