function [V, eigValEigen, eigenValueAll] = byEigen (K, d, str)
   offset = 0;

  [Q, D ] = eig((K+K')/2);

eigVal = diag(D);

[eigVal, Index] = sort(eigVal, 'descend');  Q = Q (:, Index);

if  strcmp (str, 'largest')
    
    V = Q (:, 1+offset:d+offset);
    
    eigValEigen =  eigVal (1+offset:d+offset);
    
end
if strcmp (str, 'smallest')
    
    V = Q (:, end-d+1-offset:end-offset);
    
    eigValEigen =  eigVal (end-d+1-offset:end-offset);
    
end
eigenValueAll = eigVal;
end