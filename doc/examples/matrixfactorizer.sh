MatrixFactorizer  \
  --type-space-x "lp"  \
  --px 2  \
  --type-space-y "lp"  \
  --py 2  \
  --powery 2  \
  --delta 1e-6 \
  --projection-delta 1e-6 \
  --maxiter 50 \
  --max-loops 11 \
  --orthogonal-directions 1 \
  --number-directions 1 \
  --sparse-dim 1 \
  --residual-threshold 1e-4 \
  --auxiliary-constraints "Nonnegative" \
  --max-sfp-loops 2 \
  --data pre/data.m  \
  --solution-first-factor spectralmatrix_rank2_nonnegative.m \
  --solution-second-factor pixelmatrix_rank2_nonnegative.m
