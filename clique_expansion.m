function [x_ce] = clique_expansion(incidence_list,parameter_homo_list, a, lambda, N, R)
    A = zeros(N, N);
    for i = 1:R,
        A(incidence_list{i},incidence_list{i}) = A(incidence_list{i},incidence_list{i}) + parameter_homo_list{i}/length(incidence_list{i});
    end
    for i = 1:N,
        A(i,i) = 0;
    end
    D = sum(A, 1);
    L = (1+lambda)*eye(N) - diag(D.^(-0.5))*A*diag(D.^(-0.5));
    z_ce = lambda*a.*D.^(-0.5)/L;
    x_ce = z_ce.*D.^(-0.5);
end
