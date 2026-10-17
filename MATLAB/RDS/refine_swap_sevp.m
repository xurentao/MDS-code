function [ A, Qupt, nit ] = refine_swap_sevp( A ) % nit
%REFINE_SWAP_SEVP -- for use with swap22_sevp
% A is at this stage nearly block upper triangular. We compute a
% transformation to reduce the norm of the off-diagonal block.
%
% Software based on: "Swapping 2x2 blocks in the Schur and generalized
%                     Schur form" - D. Camps, N. Mastronardi, R. Vandebril,
%                     and P. Van Dooren.
I = eye(2,2);
err = 1; % initialize to 1, only called when error is too large
Qupt = eye(4,4);
maxit = 20; nit = 0;
while (err > 10*eps) && (nit < maxit)
  nit = nit + 1;
  M = kron(I,A(3:4,3:4)) - kron(A(1:2,1:2)',I);
  b = [A(3:4,1); A(3:4,2)];
  x = M\b;
  X = reshape(x(1:4),2,2);
  [UX,SX,VX] = svd(X);SX=diag(SX);
  [COSX,SINX] = cosine_sv(SX);
  Qup = [VX*diag(COSX)*VX' VX*diag(SINX)*UX' ; ...
      -UX*diag(SINX)*VX' UX*diag(COSX)*UX'] ;
  norm_max = max(max(abs(A)));
  A = Qup'*A*Qup;
  Qupt = Qupt*Qup;
  err = max(max(abs(A(3:4,1:2))))/norm_max;
end
end

