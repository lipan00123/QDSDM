function [clustering_err, conductance, thre] = result_analysis_homo(incidence_list, parameter_homo_list, x, degree_vec, N, R)
    [x_sort, x_index] = sort(x, 'descend');
    orderinv(x_index) = 1:N;
    q = zeros(1, N);
    for i = 1:R,
        templist = orderinv(incidence_list{i});
        q(min(templist)) = q(min(templist)) + parameter_homo_list{i};
        q(max(templist)) = q(max(templist)) - parameter_homo_list{i};
    end  
    q = cumsum(q);
    conductance = inf;
    sumdegree = sum(degree_vec);
    tempnorm = 0;
    thre = 0;
    for i = 1:N-1,
        tempnorm = tempnorm + degree_vec(x_index(i));
        normalizedterm = min(tempnorm, sumdegree - tempnorm);
        if conductance > q(i)/normalizedterm,
            conductance = q(i)/normalizedterm;
            thre = i;        
        end
    end

    clustering_err = sum(x_index(1:thre) > N/2) + sum(x_index(thre+1:N) <= N/2);
    
end