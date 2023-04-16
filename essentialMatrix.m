% function ematrix = essentialMatrix (A)
%   [~, ~, V] = svd(A, 'econ');
%   v0 = V(:, end);
%   ematrix = reshape(v0, 3, 3);
%   [U, S ,V] = svd(ematrix, 'econ');
%   s = diag(S);
%   a = (s(1) + s(2))/2;
%   ematrix = U * diag([a, a, 0])* V';
%   
% end
function E = essentialMatrix (kpt1, kpt2) 
           if (size(kpt1, 2) ~= size(kpt2, 2))
                fprintf(2, 'The dimensions of the point correspondences do not match \n');
                return;
           end
            if (size(kpt1, 2) < 5)
                fprintf(2, 'Require at least 8 points for essential matrix. only %d exists\n', size(kpt1, 2));
                return;
            end
            A = normalizeImagePoints (kpt1, kpt2);
            [eigvec, ~] = byEigen (A, length(A), 'smallest');
            xopt = eigvec(:, end);
            E_opt = reshape(xopt, 3, 3);
            % rank two approxmation, no matter converged or not.
            [EU, ES, EV] = svd(E_opt);
            ES(3,3) = 0;
            a = (ES(1,1)+ ES(2,2)) / 2;
            ES(1,1) = a;
            ES(2,2) = a;
            E = EU * ES * EV';
            E =  E/norm(E, 'fro');
end
              







