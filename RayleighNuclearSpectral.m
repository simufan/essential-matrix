classdef RayleighNuclearSpectral < handle
    properties (Access = private)
        
        A = [];
        
        L = inf
        
        Lambda = 1
          
    end
  properties (Access = private)
      VectorDim = 0
        
      MatrixDimRow = 0
      MatrixDimCol = 0
        
       Lambda_L1 = 1
       Lambda_L_inf = 2
        
  end
  methods (Access = public)
      function this = RayleighNuclearSpectral (matrixA, lambda, matrix_row_dim, matrix_col_dim)
          if (nargin >= 2)
                this.A = matrixA;
                this.Lambda = lambda;
                this.Lambda_L1 = 1*lambda;
                this.Lambda_L_inf = 2.5*lambda;
            else
                fprintf(2, 'Error @RayleighNuclear object construction: The matrix A and lambda must be explicitly provided.\n');
                return;
          end
          this.VectorDim = length(matrixA);
            if (nargin == 4)
                this.MatrixDimRow = matrix_row_dim;
                this.MatrixDimCol = matrix_col_dim;
            elseif (floor(sqrt(this.VectorDim)) == sqrt(this.VectorDim))
                this.MatrixDimRow = sqrt(this.VectorDim);
                this.MatrixDimCol = sqrt(this.VectorDim);
            else
                fprintf(2, 'Error @RayleighNuclear object construction: Cannot deduce vector and matrix dimensions automatically. \n');
                fprintf(2, 'Please specific the dimension of the matrix form explictly in the constructor. \n');
                return;
            end
      end
      function [g, gradg] = gFunc (this, x)
            
            AX = this.A * x;
            
            g = x' * AX;
            
            gradg = 2 * (AX - g * x);
            
      end
      
      function [h, S] = hFunc (this, x)
            
            matx = this.vec2mat(x);
            
            [h, S] = SchattenNormRegularizer.CostL1Linf(matx, this.Lambda_L1, this.Lambda_L_inf);
            
      end
      function proxh = hFuncProx (this, x, t)
            
            matx = this.vec2mat(x);

            matx = SchattenNormRegularizer.ProximalL1Linf(matx,  t*this.Lambda_L1,  t*this.Lambda_L_inf);
            
            proxh = this.mat2vec(matx);

      end
      function L = LipschitzConst (this)
            
            this.L = 2 * svds(this.A, 1);
            
            L = this.L;
            
      end
  end
  methods (Access = private)
      function x = mat2vec (this, matx)
            if (this.VectorDim == 10 && this.MatrixDimRow == 4 && this.MatrixDimCol == 4)
                x = this.mat2vec_sym44 (matx);
            else
                x = reshape(matx, this.VectorDim, 1);
            end
      end
       function matx = vec2mat (this, x)
            if (this.VectorDim == 10 && this.MatrixDimRow == 4 && this.MatrixDimCol == 4)
                matx = this.vec2mat_sym44 (x);
            else
                matx = reshape(x, this.MatrixDimRow, this.MatrixDimCol);
            end
        end

  end
    methods  (Access = private, Static = true)
        
        function x = mat2vec_sym44 (matx)
            x = [matx(1,1); matx(1,2); matx(2,2); matx(1,3); matx(2,3); matx(3,3); matx(1,4); matx(2,4); matx(3,4); matx(4,4)];
        end
        
        
        function matx = vec2mat_sym44 (x)
            matx = [ x(1), x(2),  x(4),  x(7);
                x(2), x(3),  x(5),  x(8);
                x(4), x(5),  x(6),  x(9);
                x(7), x(8),  x(9),  x(10) ];
        end
        
    end
end