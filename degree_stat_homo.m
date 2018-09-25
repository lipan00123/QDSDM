function degree_vec = degree_stat_homo(incidence_list, parameter_homo_list, N, R)
    degree_vec = zeros(1, N);
    for i = 1:R,
        degree_vec(incidence_list{i}) = degree_vec(incidence_list{i})+(max(parameter_homo_list{i}))^2;
    end
end