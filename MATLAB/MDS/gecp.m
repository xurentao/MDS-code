function [ x, L, U, P, Q ] = gecp( A, b )
% GECP solves the linear system by Gauss elimination with complete pivoting.
% 
% Input
% -----
%   A  :  nxn matrix.
%   b  :  nx1 vector.
% Output
% ------
%   x  :  solution to Ax = b.
%   L  :  lower triangular matrix with ones as diagonals.
%   U  :  upper triangular matrix.
%   P and Q are permutations matrices so that P*A*Q = L*U 
[ n, ~ ] = size( A );
I = eye( n );
p = ( 1:n )';
q = ( 1:n )';
A_work = A;
b_work = b;
% GECP
for k = 1 : n - 1
    [ row_max, col_max ] = find( abs( A_work( k:n, k:n ) ) == max( max( abs( A_work( k:n, k:n ) ) ) ), 1 );
    m = row_max( 1 ) + k - 1; 
    l = col_max( 1 ) + k - 1; 
    % rows exchange
    if m ~= k
        A_work( [ k, m ], : ) = A_work( [ m, k ], : );
        b_work( [ k, m ]) = b_work( [m, k] );
        p( [ k, m ] ) = p( [ m, k ] );
    end
    % columns exchange
    if l ~= k
        A_work( :, [ k, l ] ) = A_work( :, [ l, k ] );
        q( [ k, l ] ) = q( [ l, k ] );
    end
    A_work( k+1:n, k ) = A_work( k+1:n, k ) / A_work( k, k );
    A_work( k+1:n, k+1:n ) = A_work( k+1:n, k+1:n ) - A_work( k+1:n, k ) * A_work( k, k+1:n );
end
L = eye(n) + tril( A_work, -1 );
U = triu( A_work );    
% solve Ly = P*b
y = zeros( n, 1 );
for i = 1:n
    y( i ) = b_work( i ) - L( i, 1:i-1 ) * y( 1:i-1 );
end
% solve U*z = y
z = zeros( n, 1 );
for i = n:-1:1
    z( i ) = ( y( i ) - U( i, i+1:n ) * z( i+1:n ) ) / U( i, i );
end 
% obtain x = Q*z
x = zeros( n, 1 );
x( q ) = z;
P = I( p, : );
Q = I( :, q );
end