classdef TwoViewGeometry
    methods(Access = public, Static = true)
         function [E, x, y1, y2] = EMatrix (kpt1, kpt2)

                ref_lambda = 0.1;
                [E, ~, x, y1, y2] = TwoViewGeometry.EMatrixPGS (kpt1, kpt2, 10, ref_lambda);
                 %E = TwoViewGeometry.calculateEMatrix (kpt1, kpt2, ref_lambda);
 
         end
     
         function [E, info_pgs,x, y1, y2] = EMatrixPGS (kpt1, kpt2, iters, lambda)
                options.debug = false;
                options.max_iters = iters;
                options.method = 'PGS'; % AMPGS, APGS
                options.error_tol = 1e-5;
                if ~exist('lambda', 'var') || (exist('lambda', 'var') && isempty(lambda))
                    ref_lambda = 0.005;
                else
                    ref_lambda = lambda;
                end
                [E, info_pgs, x, y1, y2] = TwoViewGeometry.calculateEMatrix (kpt1, kpt2, ref_lambda, options);
         end

         function [E, info_pgs] = EMatrixPGSDebug (kpt1, kpt2, lambda, options)
                if ~exist('lambda', 'var') || (exist('lambda', 'var') && isempty(lambda))
                    ref_lambda = 0.005;
                else
                    ref_lambda = lambda;
                end
                if ~exist('options', 'var') || (exist('options', 'var') && ~isfield(options, 'debug'))
                    options.debug = true;
                end
                if ~isfield(options, 'method')
                    options.method = 'PGS'; % PGS, AMPGS, APGS
                end
                if ~isfield(options, 'max_iters')
                    options.max_iters = 50;
                end
                if ~isfield(options, 'error_tol')
                    options.error_tol = 1e-8;
                end
                [E, info_pgs] = TwoViewGeometry.calculateEMatrix (kpt1, kpt2, ref_lambda, options);
         end
          
          
        
         function pts3D = TrangulatePoint3D (P1, P2, kpt1, kpt2)
            
            % this is an opensource implementation from Manolis Lourakis
            pts3D = stereoReconsPts(P1, P2, kpt1, kpt2);
            
            % this requires computer vision toolbox
            % pts3D=triangulate(kpt1', kpt2', P1', P2')';

         end
        
    end
   methods(Access = private, Static = true)
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
     
       function [E, info, x, y1, y2] = calculateEMatrix (kpt1, kpt2, lambda, options) 
           if (size(kpt1, 2) ~= size(kpt2, 2))
                fprintf(2, 'The dimensions of the point correspondences do not match \n');
                return;
           end
            if (size(kpt1, 2) < 5) 
                fprintf(2, 'Require at least 8 points for essential matrix. only %d exists\n', size(kpt1, 2));
                return;
            end
            [A] = normalizeImagePoints (kpt1, kpt2);
            [eigvec, eigval] = TwoViewGeometry.byEigen (A, length(A), 'smallest');
            x0 = eigvec(:, end);
            LipschtzConst = 2 * max(abs(eigval));
            if (lambda == 0)
                xopt = x0;
            else
                pgs_problem = RayleighNuclearSpectral(A, lambda);
                pgs_method = ProxGradSphere (pgs_problem);
                pgs_method.MaxIterations = options.max_iters;
                pgs_method.ErrorThreshold = options.error_tol;
                if (options.debug)
                    pgs_method.DebugFlag = true;
                end
                if strcmp(options.method, 'PGS')
                    if options.debug
                        [xopt, info] = pgs_method.Optimize (x0,  1.0/LipschtzConst);
                    else
                        [xopt, info, x, y1, y2] = pgs_method.OptimizeByNesteroveMomentum(x0,  1.0/LipschtzConst, false, false);
                    end
                end
                if strcmp(options.method, 'AMPGS')
                    [xopt, info, x, y1, y2] = pgs_method.OptimizeByNesteroveMomentum(x0,  1.0/LipschtzConst, true, true);
                end
                if strcmp(options.method, 'APGS')
                    [xopt, info, x, y1, y2] = pgs_method.OptimizeByNesteroveMomentum(x0,  1.0/LipschtzConst, true, false);
                end
            end
            E_opt = reshape(xopt, 3, 3);
            % rank two approxmation, no matter converged or not.
            [EU, ES, EV] = svd(E_opt);
            ES(3,3) = 0;
            a = (ES(1,1)+ ES(2,2)) / 2;
            ES(1,1) = a;
            ES(2,2) = a;
            E = EU * ES * EV';
            E =  E/norm(E, 'fro');
            info.LipschtzConst = LipschtzConst;
                   
       end
    
   end

 
end   
function matr3by3 = skewSym (v)
            
            matr3by3 = [  0,     -v(3),   v(2);
                                    v(3),    0,     -v(1);
                                   -v(2),   v(1),    0    ];
            
end