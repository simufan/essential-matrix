classdef ProxGradSphere < handle
    properties (Access = public)
        MaxProxyStepSize = inf
        
        AdaptiveMaxProxyStepSize = false;  % if true. Each iteration, we set MaxProxyStepSize to the previous working proxyStepsize.
        ShrinkingFactor = 0.8
        
        MaxIterations = 50
        Verbose = true
        
        DebugFlag = false
        ErrorThreshold = 1e-3
        
        CostFunctionG = []        
        CostFunctionH = []
        ConvergenceDeltaXopt = []
        ConvergenceProjTx_U = []
        ConvergenceTangentV = []
        ConvergenceDeltaFcost = []
    end
    
    properties (Access = private)
        
        LipschitzConst = inf
    
        ProblemHandle
        
    end
    methods  (Access = public)
        
        function this = ProxGradSphere (ProblemInstance)
        
            this.ProblemHandle = ProblemInstance;
            
        end
        function [x_opt, info,x, y1, y2] =  OptimizeByNesteroveMomentum (this, xInit, maxProxyStepSize, useNesterovMomentumAcceleration, useBeckerMonotoneVersion)
            xInit = xInit ./ norm(xInit, 'fro');

            if ~exist('maxProxyStepSize', 'var') || (exist('maxProxyStepSize', 'var') && isempty(maxProxyStepSize))
                maxProxyStepSize = -1;
            end
            if ~exist('useNesterovMomentumAcceleration', 'var') || (exist('useNesterovMomentumAcceleration', 'var') && isempty(useNesterovMomentumAcceleration))
                useNesterovMomentumAcceleration = true;
            end
            if ~exist('useBeckerMonotoneVersion', 'var') || (exist('useBeckerMonotoneVersion', 'var') && isempty(useBeckerMonotoneVersion))
                useBeckerMonotoneVersion = true;
            end
            this.CostFunctionG = [];
            this.CostFunctionH = [];

            this.ConvergenceDeltaXopt = [];
            this.ConvergenceProjTx_U = [];
            this.ConvergenceTangentV = [];
            this.ConvergenceDeltaFcost = [];
            initProxyStepsizeIter = 0;

            % if maxStepSize  is non-positive (or infinity), we find a working maxStepSize by adaptive line-search
            if (maxProxyStepSize <= 0 || maxProxyStepSize == inf)
                [maxProxyStepSize, initProxyStepsizeIter] = this.AdaptiveSearchFirstProxyStepSize (xInit, inf);
            end
            this.MaxProxyStepSize = maxProxyStepSize;
            usedIter = 0;
            
            totalIter = 0;
            
            line_search_cnt = 0;
            auxiliary_weight = 1;
            
            auxiliary_state = xInit;
            
            x_prev  = xInit;
            g_prev = this.ProblemHandle.gFunc(x_prev);
            [h_prev, S_prev] =  this.ProblemHandle.hFunc(x_prev);
            norm_projTx_U = nan;
            proxyStepsize = min(this.MaxProxyStepSize, 0.7/this.ProblemHandle.hFunc(auxiliary_state));
            this.CostFunctionG = [this.CostFunctionG, g_prev];
            this.CostFunctionH = [this.CostFunctionH, h_prev];
            if (this.Verbose)
%                 fprintf(1, '[%d]   g+h = %f,  g = %f,  h = %f, |projTx_U| = %f | LSIter = %d | proxstepsize = %f ', usedIter, g_prev+h_prev, g_prev, h_prev, norm_projTx_U,  line_search_cnt, proxyStepsize);
                if (1)
%                     fprintf(1, '\n S = ['); fprintf(1, '%f ', S_prev); fprintf(1, ']');
                end
                fprintf(1, '\n');
            end
            %plot
             x = zeros(this.MaxIterations, 0);
             y1 = zeros(this.MaxIterations, 0);
             y2 = zeros(this.MaxIterations, 0);
            while usedIter < this.MaxIterations
                line_search_cnt =  line_search_cnt + 1;
                if ~useNesterovMomentumAcceleration
                    auxiliary_weight = 1;
                end
                totalIter = totalIter + 1;
                [g, gradg] = this.ProblemHandle.gFunc(auxiliary_state);
                
                proxh = this.ProblemHandle.hFuncProx(auxiliary_state - proxyStepsize * gradg, proxyStepsize);

                cconst = 1 / (auxiliary_state' * proxh);
                t = cconst * proxyStepsize;  % real step size: t
                
                v = cconst * proxh - auxiliary_state; % tangent update: v
                projTx_U = v/t;
                norm_projTx_U = norm(projTx_U);
                norm_tangent_V = norm(v);
                x_opt = (cconst / norm(cconst * proxh)) * proxh;
                g_opt = this.ProblemHandle.gFunc(x_opt);
                norm_deltaXopt = norm(x_opt - x_prev, 'fro');
                rhs =  g +gradg' * v + (v' * v)/(2*t);
                % line-search condition to verify for the next estimate
                lineSearchCond = (g_opt <=  rhs );
                if lineSearchCond
                    usedIter = usedIter + 1;
                    auxiliary_weight_next = (sqrt(4*auxiliary_weight*auxiliary_weight + 1) + 1)/2 ; % this is accelerated version
                    g_opt = this.ProblemHandle.gFunc(x_opt);
                    [h_opt, S_opt] = this.ProblemHandle.hFunc(x_opt);
                    f_descent =  (g_prev + h_prev) - (g_opt + h_opt);
                    if useNesterovMomentumAcceleration && useBeckerMonotoneVersion && (f_descent < 0)
                        % difference vector
                        dv = (auxiliary_weight/auxiliary_weight_next) * this.invRetr(x_prev, x_opt);
                        x_opt = x_prev;
                        g_opt = g_prev;
                        h_opt = h_prev;
                        S_opt = S_prev;
%                         fprintf(2, 'a non-monotone step\n');
                    else
                        % difference vector
                        dv = ((1-auxiliary_weight)/auxiliary_weight_next) * this.invRetr(x_opt, x_prev);
                    end
                    auxiliary_state =  this.Retr(x_opt,  dv);
                    auxiliary_weight = auxiliary_weight_next;     % set this one to 1, we obtain standard PGS algorithm.  auxiliary_state = x_prev
                     x(usedIter) = usedIter;
                     y1(usedIter) = S_opt(1);
                     y2(usedIter) = S_opt(2);
                    if (this.Verbose)
%                         fprintf(1, '[%d]   g+h = %f,  g = %f,  h = %f,  f_descent = f_{k} - f_{k+1} = %.20f,  ||projTx_U|| = %f,  ||deltaXopt|| = %f,  ||tangentV|| = %f,  LSIter = %d,  proxstepsize = %f,  stepsize = %f ', usedIter, g_opt+h_opt, g_opt, h_opt, f_descent, norm_projTx_U, norm_deltaXopt,  norm_tangent_V, line_search_cnt, proxyStepsize, t);
                        if (1)
%                             fprintf(1, '\n S = ['); fprintf(1, '%f ', S_opt); fprintf(1, ']');
                          
                        end
                        fprintf(1, '\n');
                    end
                    line_search_cnt = 0;
                    
                     this.CostFunctionG = [this.CostFunctionG, g_opt];
                     this.CostFunctionH = [this.CostFunctionH, h_opt];

                     this.ConvergenceDeltaXopt = [this.ConvergenceDeltaXopt,  norm_deltaXopt];
                     this.ConvergenceProjTx_U = [this.ConvergenceProjTx_U,  norm_projTx_U];
                     this.ConvergenceTangentV = [this.ConvergenceTangentV,  norm_tangent_V];
                     this.ConvergenceDeltaFcost = [this.ConvergenceDeltaFcost,  f_descent];
                     x_prev = x_opt;
                    g_prev = g_opt;
                    h_prev = h_opt;
                    S_prev = S_opt;
                    if (this.AdaptiveMaxProxyStepSize)
                        this.MaxProxyStepSize = proxyStepsize;
                    end
                                        
                    proxyStepsize = min(this.MaxProxyStepSize, 0.7/this.ProblemHandle.hFunc(auxiliary_state));
                    if (abs(f_descent) < 1e-18) || ( norm_projTx_U < this.ErrorThreshold)
                        break;
                    end
               else
 
               proxyStepsize = proxyStepsize * this.ShrinkingFactor;
               end
            end
%             fprintf(1, 'usedIter = %d,  lineSearchIter = %d, totalIter = %d \n\n', usedIter, totalIter-usedIter, totalIter);
            
            info.initProxyStepsizeIter = initProxyStepsizeIter;

            info.totalIter = totalIter;
            info.usedIter = usedIter;
            info.CostFunctionArr = [ this.CostFunctionG+this.CostFunctionH;
                                                        this.CostFunctionG; 
                                                        this.CostFunctionH; ];
            info.ConvergenceDeltaXopt = this.ConvergenceDeltaXopt;
            info.ConvergenceProjTx_U = this.ConvergenceProjTx_U;
            info.ConvergenceTangentV = this.ConvergenceTangentV;
            info.ConvergenceDeltaFcost = this.ConvergenceDeltaFcost;
        end
        function [proxyStepsize, iters] = AdaptiveSearchFirstProxyStepSize (this, x_cur, maxProxyStepsize)
            if ~exist('maxProxyStepsize', 'var') || (exist('maxProxyStepsize', 'var') && isempty(maxProxyStepsize))
                maxProxyStepsize = inf;
            end
            x_cur = x_cur ./ norm(x_cur, 'fro');

            ub = 0.7 / this.ProblemHandle.hFunc(x_cur);
            proxyStepsize = min(abs(maxProxyStepsize),  ub);

            [g, gradg] = this.ProblemHandle.gFunc(x_cur); 

            iters = 0;
            found = false;
            while (iters < 1000)
                iters = iters + 1;

                proxh = this.ProblemHandle.hFuncProx(x_cur - proxyStepsize * gradg, proxyStepsize);

                cconst = 1 / (x_cur' * proxh);
                t = cconst * proxyStepsize;

                v = cconst * proxh - x_cur;
                x_next = (cconst / norm(cconst * proxh)) * proxh;

                lineSearchCond = ( this.ProblemHandle.gFunc(x_next) <= g + gradg' * v + (v' * v)/(2*t) );
                if (lineSearchCond)
                    if (proxyStepsize == ub)
                        return;
                    end
                    found = true;
                    proxyStepsize = min(2 * proxyStepsize,  ub);
                else
                    if (found)
                        proxyStepsize = 0.5 * proxyStepsize;
                        return;
                    end
                    proxyStepsize = 0.1 * proxyStepsize;

                end
            end
        end
        
        function [x_opt, info] = Optimize (this, xInit, maxProxyStepSize)
            if ~exist('maxProxyStepSize', 'var') || (exist('maxProxyStepSize', 'var') && isempty(maxProxyStepSize))
                maxProxyStepSize = -1;
            end
            xInit = xInit ./ norm(xInit, 'fro');
            SingularValuesArr = [];
            CostFunctionArr = [];
            if (maxProxyStepSize <= 0 || maxProxyStepSize == inf)
                maxProxyStepSize = this.AdaptiveSearchFirstProxyStepSize (xInit, inf);
            end
            this.MaxProxyStepSize = maxProxyStepSize;
            x_cur = xInit;
            
            usedIter = 0;
            
            totalIter = 0;
            IntervalWidth = 1/ this.ProblemHandle.hFunc(x_cur);
            proxyStepsize = min(this.MaxProxyStepSize, IntervalWidth);
            norm_projTx_U = nan;
          
            while usedIter < this.MaxIterations
                totalIter = totalIter + 1;
                
                [g, gradg] = this.ProblemHandle.gFunc(x_cur);
                
                proxh = this.ProblemHandle.hFuncProx(x_cur - proxyStepsize * gradg, proxyStepsize);
                cconst = 1 / (x_cur' * proxh);
                
                t = cconst * proxyStepsize;
                v = cconst * proxh - x_cur;
                
                projTx_U = v/t;
                norm_projTx_U = norm(projTx_U);
                x_opt = (cconst / norm(cconst * proxh)) * proxh;
                
                lineSearchCond = ( this.ProblemHandle.gFunc(x_opt) <= g + gradg' * v + (v' * v)/(2*t) );
                IntervalWidth = 1/ this.ProblemHandle.hFunc(x_cur);
                if lineSearchCond
                    if(this.DebugFlag && usedIter < 5)
                        test_proxy_t =  1.0*(-0.97 : 0.01 : 0.97) * IntervalWidth;
                        test_t = zeros(1, length(test_proxy_t));
                        test_c = test_t; 
                        for ii = 1 : length(test_proxy_t)
                            y_star = this.ProblemHandle.hFuncProx(x_cur - test_proxy_t(ii) * gradg, abs(test_proxy_t(ii)));
                            test_c(ii) = x_cur' * y_star;
                            test_t(ii) = test_proxy_t(ii) / test_c(ii);
                        end
                        info.StepsizeMapping{usedIter+1} = [test_proxy_t; test_t];
                        info.StepsizeMapping_c{usedIter+1} = [test_proxy_t; test_c];
                    end
                    if (this.Verbose)
                        [h, S] = this.ProblemHandle.hFunc(x_cur);
%                         fprintf(1, '[%d]   g+h = %f,  g = %f,  h = %f, |projTx_U| = %f,  S = [', usedIter, g+h, g, h, norm_projTx_U); fprintf(1, '%f ', S); fprintf(1, ']\t');
%                         fprintf(1, 'proxyStepsize = %f \n', proxyStepsize);
                    end
                    [h, S] = this.ProblemHandle.hFunc(x_cur);
                    SingularValuesArr = [ SingularValuesArr,  [S(1); S(2); S(3)] ];
                    CostFunctionArr = [ CostFunctionArr,  [g+h; g; h] ];
                    x_cur = x_opt;
                    
                    usedIter = usedIter + 1;
                    if (this.AdaptiveMaxProxyStepSize)
                        this.MaxProxyStepSize = proxyStepsize;
                    end
                    proxyStepsize = min(this.MaxProxyStepSize, IntervalWidth);
                    if ( norm_projTx_U < this.ErrorThreshold)  % || ( abs(g - g_opt) < 1e-3 )
                        break;
                    end
                else
                    proxyStepsize = proxyStepsize * this.ShrinkingFactor;
 
                end
            end
            g_opt = this.ProblemHandle.gFunc(x_opt);
            [h_opt, S_opt] = this.ProblemHandle.hFunc(x_opt);
            SingularValuesArr = [ SingularValuesArr,  [S_opt(1); S_opt(2); S_opt(3)] ];
            CostFunctionArr = [ CostFunctionArr,  [g_opt+h_opt; g_opt; h_opt] ];
            info.usedIter = usedIter;
            info.CostFunctionArr = CostFunctionArr;
            info.SingularValuesArr = SingularValuesArr;
            if (this.Verbose)
                fprintf(1, '[%d]   g+h = %f,  g = %f,  h = %f,  S = [', usedIter, g_opt+h_opt, g_opt, h_opt); fprintf(1, '%f ', S_opt); fprintf(1, ']\n\n');
                fprintf(1, 'usedIter = %d,  totalIter = %d \n\n', usedIter, totalIter);
            end
        end
    end
    
    methods (Access = private, Static = true)
        
        % the v that satisfies:  x_next = R_x (v)
        function v = invRetr (x_cur, x_next)
            c = 1 / (x_cur' * x_next);
            v = c * x_next - x_cur;
        end
        
        function x_next = Retr(x_cur, v)
            x_next = (x_cur + v);
            x_next = x_next / norm(x_next);
        end
        
    end
end
