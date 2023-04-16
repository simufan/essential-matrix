classdef SchattenNormRegularizer
    methods (Access = public, Static = true)
        function [h, S] = CostL1Linf (matx, lambda_1, lambda_inf)
            
            S = svd(matx);
            
            h = lambda_1 * sum(S) + lambda_inf * S(1);
            
        end
        function proxh = ProximalL1Linf  (matx, t_1,  t_inf)
            [U, S, V] = svd (matx, 'econ');
            
            vs = diag(S);
            
            vs = vs - t_1;  % Schattern L1 norm, nuclear norm thresholding
            vs(1) = vs(1) - t_inf; % Schattern L-infity norm, spectral norm thresholding
            
            vs = max(vs, 0);
            
            proxh = U * diag(vs) * V';
                     
   
            
        end
    end
    


end