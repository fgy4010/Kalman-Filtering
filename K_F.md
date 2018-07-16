# Kalman-Filtering
Java to implement Kalman Filtering

    public int predictValueKalman(int p, int q, int period)
    {
      double [] data = this.preDealDiff(period);
      int m = Math.max(p, q + 1);
      // calculate T
      double [][] T = new double[m][m];       //default initialization.
      if (p == 0)
      {
          for (int i = 0; i < m; ++i)
          {
              T[i][0] = 0;
          }
      }
      else
      {
          double [] ARCoe = this.arimaCoe.get(0);
          for (int i = 0; i < m; ++i)
          {
              if (i < ARCoe.length - 1)
              {
                  T[i][0] = ARCoe[i];
              }
              else
              {
                  T[i][0] = 0;
              }
          }
      }
      for (int i = 0; i < m; ++i)
      {
          if (i + 1 < m)
          {
              T[i][i + 1] = 1;
          }
          else
          {
              T[i][0] = 0;
          }
      }
      Matrix TMarix = new Matrix(T);      //calcilate matrix    
      
      // calculate R
      double [] R = new double[m];       
      double [] theta = new double[m - 1];
      R[0] = 1;
      if (q == 0)
      {
          for (int i = 1; i < m; ++i)
          {
              R[i] = 0;
          }
      }
      else
      {
          double [] MACoe = null;
          if (p == 0)
          {
              MACoe = this.arimaCoe.get(0);
          }
          else
          {
              MACoe = this.arimaCoe.get(1);
          }
          for (int i = 1; i < m; ++i)
          {
              if (i < MACoe.length)
              {
                  R[i] = MACoe[i];
                  theta[i - 1] = MACoe[i];
              }
              else
              {
                  R[i] = 0;
              }
          }
      }
      Matrix RMatrix = new Matrix(R, 1).transpose();
      //calculate Z
      double [] Z = new double[m];        
      Z[0] = 1;
      Matrix ZMatrix = new Matrix(Z, 1);

      //calculate Q
      double Q = 0.0;
      if (p == 0)     //MA
      {
          Q = this.arimaCoe.get(0)[0];
      }
      else if (q == 0)        //AR
      {
          double [] ARCoe = this.arimaCoe.get(0);
          Q = ARCoe[ARCoe.length - 1];
      }
      else        //ARMA
      {
          Q = this.arimaCoe.get(1)[0];
      }
      Q = Math.sqrt(Q);
      
      //calculate Matrix sigma
      double [] noise = new double[data.length];
      for (int i = 0; i < noise.length; ++i)
      {
          noise[i] = new Random().nextGaussian() * Q;
      }
      int len = data.length - m;
      double [][] simData = new double[len][m];
      for (int i = 0; i < len; ++i)
      {
          for (int j = 0; j < m; ++j)
          {
              double ar = 0.0;
              for (int k = j; k < m; ++k)
              {
                  ar += T[k][0] * data[(m - 1 - k + i) + j];
              }
              double ma = 0.0;
              for (int k = j; k < m - 1; ++k)
              {
                  ma += theta[k] * noise[(m - 1 - k + i) + j];
              }
              simData[i][j] = ar + ma;
          }
      }
      double [][] sigmaTrans = new double[m][len];
      for (int i = 0; i < m; ++i)
      {
          for (int j = 0; j < len; ++j)
          {
              sigmaTrans[i][j] = simData[j][i];
          }
      }
      double [][] sigma = new ARMAMaths().computeMutalCorrMatrix(sigmaTrans);
      Matrix sigmaMatrix = new Matrix(sigma).transpose();
      
      // Kalman iterator
      // end loop
      double [] v = new double[data.length];
      double [][] xhat = new double[m][data.length];
      double [] tmphat = new double[m];
      double V = 0.0;
      double [] K = new double[m];
      double [][] L = new double[m][m];
      Matrix xhatMatrix =  new Matrix(tmphat, 1).transpose();
      Matrix KMatrix = new Matrix(K, 1).transpose();
      Matrix LMatrix = new Matrix(L);
      for (int i = 1; i < data.length; ++i)
      {
          Matrix tmp = ZMatrix.times(xhatMatrix);
          v[i] = data[i - 1] - tmp.get(0, 0);
          V = ZMatrix.times(sigmaMatrix).times(ZMatrix.transpose()).get(0, 0);
          KMatrix = TMarix.times(sigmaMatrix).times(ZMatrix.transpose()).times(1.0 / V);
          LMatrix = TMarix.minus(KMatrix.times(ZMatrix));
          xhatMatrix = TMarix.times(xhatMatrix).plus(KMatrix.times(v[i]));
          for (int k = 0; k < xhatMatrix.getRowDimension(); ++k)
          {
              xhat[k][i] = xhatMatrix.get(k, 0);
          }
          Matrix tmpMatrix = RMatrix.times(Q).times(RMatrix.transpose());
          sigmaMatrix = TMarix.times(sigmaMatrix).times(LMatrix.transpose()).plus(tmpMatrix);
      }
      int predictValue = (int)Math.round(xhat[0][xhat[0].length - 1]);
      return predictValue;
  }
