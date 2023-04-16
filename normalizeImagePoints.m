   function [A] = normalizeImagePoints(kpt1, kpt2)
           n = size(kpt1, 2);
           k = [2759.48 0 1520.69; 
                0 2764.16 1006.81;
                0 0 1];
            nkpt1 = inv(k) * kpt1;

            nkpt2 = inv(k) * kpt2;

            A = zeros(9, 9);
            for ii = 1 : n
                p1 = [nkpt1(1:2,ii); 1];
                p2 = [nkpt2(1:2,ii); 1];
                v = kron(p1', p2');
                A = A + v' * v;
            end

            A = A / n;
       end

