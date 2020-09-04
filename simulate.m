
%% Simulation settings
Pool_mat=kirkman_42_238;
N=238; M=42;
num_populations=2;
n_sick=6;

%% Generate populations
xs=zeros(N,num_populations);
for n = 1:num_populations
    positive_samples=zeros(N,1);
    positive_samples(1:n_sick)=1;
    positive_samples=positive_samples(randperm(N));
    x_true=1000*rand(N,1).*(positive_samples);
    xs(:,n)=x_true;
end

%% Simulate Poisson pooled results
y_means=zeros(M,num_populations);
ys=zeros(M,num_populations);
for n=1:num_populations,
    x_true=xs(:,n);
    Pool_means=Pool_mat*x_true;
    y_means(:,n)=Pool_means;
    Pool_values=poissrnd(Pool_means);
    ys(:,n)=Pool_values;
end
sum(sum(y_means>0))/num_populations;
%% CS with entropy reg
xs_estimate=zeros(N,num_populations);
alpha=30;
for n = 1:num_populations
    Pool_values=ys(:,n);
    pool_distr=ones(M,N);
    pool_distr(min(bsxfun(@times,Pool_mat,Pool_values)<1E-4,(Pool_mat>0)))=0;
    support=find(prod(pool_distr,1));
    %support=1:N %Uncomment if optimized_mat is used
    num_iter=4;
    x_estimate=decode_entropy(Pool_values,Pool_mat,alpha,num_iter,support,N);
    xs_estimate(:,n)=x_estimate;
end
xs_entropy=xs_estimate;
figure()
scatter(xs(:),xs_entropy(:));
hold;
plot([0,max(xs(:))],[0,max(xs(:))]);
xlabel('True viral concentration');
ylabel('Estimated viral concentrations');


%% CS with l1/2 reg
xs_estimate=zeros(N,num_populations);
alpha=5;
for n = 1:num_populations
    Pool_values=ys(:,n);
    pool_distr=ones(M,N);
    pool_distr(min(bsxfun(@times,Pool_mat,Pool_values)<1E-4,(Pool_mat>0)))=0;
    support=find(prod(pool_distr,1));
    %support=1:N %Uncomment if optimized_mat is used
    num_iter=8;
    x_estimate=decode_sparse(Pool_values,Pool_mat,alpha,num_iter,support,N);
    xs_estimate(:,n)=x_estimate;
end
figure()
xs_sparse=xs_estimate;
scatter(xs(:),xs_sparse(:))
hold;
plot([0,max(xs(:))],[0,max(xs(:))]);
xlabel('True viral concentration')
ylabel('Estimated viral concentrations')
title('l1/2 penalty')


%% cs with l1reg
xs_estimatel1=zeros(N,num_populations);
alpha=10;
for n = 1:num_populations
    Pool_values=ys(:,n);
    pool_distr=ones(M,N);
    pool_distr(min(bsxfun(@times,Pool_mat,Pool_values)<1E-4,(Pool_mat>0)))=0;
    support=find(prod(pool_distr,1));
    %support=1:N %Uncomment if optimized_mat is used
    x_estimate=decode_l1(Pool_values,Pool_mat,alpha,support,N);
    xs_estimatel1(:,n)=x_estimate;
end
figure()
scatter(xs(:),xs_estimatel1(:))
hold;
title('l1 penalty')
xlabel('True viral concentration')
ylabel('Estimated viral concentrations')
x_l1=xs_estimatel1;