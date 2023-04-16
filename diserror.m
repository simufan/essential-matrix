        function [error,in] = diserror(E, kpt1, kpt2)
                n = size(kpt1, 2);
                in = zeros(n, 1);
                rsd = 0;
                for ii = 1 : n
                    x1 = kpt1(:,ii);
                    x2 = kpt2(:,ii);
                    l = E * x1 ;
                    a = x2' * l;
                    rsd = rsd + sqrt(a * a / (l(1) * l(1) + l(2) * l(2))) ;
                end
                 error = rsd/n;
        end