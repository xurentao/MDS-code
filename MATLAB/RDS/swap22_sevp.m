function [A,Q,tol] = swap22_sevp(A, qrup)
% Swap 2x2 blocks in real Schur form.
% 
% Input
% -----
%   A  :  4x4 block upper triangular having 2x2 blocks with complex-conjugate
%         eigenvalues.
%   qrup :  boolean, if true, the transformations are computed via a QR
%         decomposition, otherwise via an SVD.
%
% Output
% ------
%   A  :  4x4 block upper triangular with swapped eigenvalues
%   Q  :  similarity transformation
%   tol  :  error indicator 
%
% Software based on: "Swapping 2x2 blocks in the Schur and generalized
%                     Schur form" - D. Camps, N. Mastronardi, R. Vandebril,
%                     and P. Van Dooren.
warning('off','MATLAB:nearlySingularMatrix');
I = eye(2,2);
M = kron(I,A(1:2,1:2)) - kron(A(3:4,3:4)',I);
x = M\[A(1:2,3); A(1:2,4)];
X = reshape(x(1:4),2,2);
if qrup % Use QR decomposition
[Q,~] = qr([-X; I]);
else % Use SVD
[UX,SX,VX] = svd(X);SX=diag(SX);
[COSX,SINX] = cosine_sv(SX);
Q = [-UX*diag(SINX)*VX' UX*diag(COSX)*UX' ; ...
     VX*diag(COSX)*VX' VX*diag(SINX)*UX'] ;   
end
norm_max = max(max(abs(A)));
A = Q'*A*Q;
tol = max(max(abs(A(3:4,1:2))))/norm_max;
if tol > 10*eps
  [A,Qup] = refine_swap_sevp(A);
  Q = Q * Qup;
end
A(3:4,1:2) = 0; 
warning('on','MATLAB:nearlySingularMatrix');
end