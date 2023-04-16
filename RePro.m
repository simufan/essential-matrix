function reprojError = RePro(FMat, kpt1, kpt2)
    FMat = FMat / norm(FMat, 'fro');
    P1 = [ eye(3),   zeros(3, 1) ];
    epole2 = nullVector(FMat');
    P2 = [skewSym(epole2) * FMat,  epole2];
    pts3D=triangulate(kpt1', kpt2', P1', P2')';
    reprojError = ReError(kpt1, kpt2, P1, P2, pts3D);
end

       function v = nullVector (E)

                [~,~,V] = svd (E);

                v = V(:,end);

        end

        function e = ReError(q1_j, q2_j, P1, P2, pts3D)
                m = length(pts3D);
                Q_j = [pts3D; ones(1, m) ];
                E = ReErrorRes(q1_j, q2_j, P1, P2, Q_j);
                e = sqrt(sum(sum(E.^2)) / (m*2));
        end

        function [E, ohatq1_j, ohatq2_j] = ReErrorRes(q1_j, q2_j, P1, P2, Q_j)
                hatq1_j = P1 * Q_j;
                hatq1_j = [ hatq1_j(1, :)./hatq1_j(3, :) ; hatq1_j(2, :)./hatq1_j(3, :) ];
                hatq2_j = P2 * Q_j;
                hatq2_j = [ hatq2_j(1, :)./hatq2_j(3, :) ; hatq2_j(2, :)./hatq2_j(3, :) ];
                E = [
                    q1_j - hatq1_j
                    q2_j - hatq2_j
                    ];
                if nargout > 1, ohatq1_j = hatq1_j; end
                if nargout > 2, ohatq2_j = hatq2_j; end
        end
        
        function matr3by3 = skewSym (v)

                matr3by3 = [           0,        -v(3),   v(2);
                                                   v(3),         0,           -v(1);
                                                  -v(2),   v(1),          0         ];

          end


