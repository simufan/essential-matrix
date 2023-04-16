 function error = reprojection(kpt1, kpt2, K, pts3d, R, t)
           n = size(kpt1, 2);
           pts3d = [pts3d; ones(1, n)];
           P1 = K * [ eye(3),  zeros(3, 1)];
           P2 = K * [ R, t];
           error = 0;
           for i = 1 : n
                p1 = P1 * pts3d(:, i);
               
                p2 = P2 * pts3d(:, i);
               
                p1 =  [p1(1)/p1(3) ; p1(2)/p1(3)];
                
                p2=   [p2(1)/p2(3) ; p2(2)/p2(3)];
                
                error = error +(p1(1)- kpt1(1,i))^2 + (p1(2)- kpt1(2,i))^2 + (p2(1)- kpt2(1,i))^2+ (p2(2)- kpt2(2,i))^2;
           end
           error = sqrt(error/(2 * n));
        end