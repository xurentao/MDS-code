function [ A_work, Q, tol ] = mds( A ) 
% Swap 2-by-2 blocks in the real Schur form by MDS.
% 
% Input
% -----
%   A  :  4-by-4 block upper triangular having 2-by-2 blocks with complex-conjugate
%         eigenvalues.
%   wanted  :  boolean, if ture, output transformation Q, otherwise Q is the 4-by-4 identity matrix.
% Output
% ------
%   A_work  :  4-by-4 block upper triangular with swapped eigenvalues.
%   tol  :  error metric
%   Q  :  orthogonal similarity transformation.
A_work = A; I = eye( 2 ); Q = eye( 4 );
A11 = A( 1:2, 1:2 ); A22 = A( 3:4, 3:4 );
A12 = A( 1:2, 3:4 ); b = reshape( A12, [ 4, 1 ] );
T = kron( I, A11 ) - kron( A22', I );
[ x, ~, ~, ~, ~ ] = gecp( T, b );
X = reshape( x, [ 2, 2 ] );
[ ~, col ] = find( abs( X( 1:2, 1:2 ) ) == max( max( abs( X( 1:2, 1:2 ) ) ) ), 1 );
H = [ -X; I ];
% pivoting Givens QR
if col == 1
    [ G1, H( [ 1, 2 ], 1 ) ] = planerot( H( [ 1, 2 ], 1 ) ); 
    H( [ 1, 2 ], 2 ) = G1 * H( [ 1, 2 ], 2 ); 
    [ G2, ~ ] = planerot( H( [ 1, 3 ], 1 ) );
    H( [ 1, 3 ], 2 ) = G2 * H( [ 1, 3 ], 2 );
    [ G3, H( [ 2, 4 ], 2 ) ] = planerot( H( [ 2, 4 ], 2 ) );
    [ G4, ~ ] = planerot( H( [ 2, 3 ], 2 ) );

    A_work( [ 1, 2 ], : ) = G1 * A_work( [ 1, 2 ], : );
    A_work( :, [ 1, 2 ] ) = A_work( :, [ 1, 2 ] ) * G1';
    A_work( [ 1, 3 ], : ) = G2 * A_work( [ 1, 3 ], : );
    A_work( :, [ 1, 3 ] ) = A_work( :, [ 1, 3 ] ) * G2';
    A_work( [ 2, 4 ], : ) = G3 * A_work( [ 2, 4 ], : );
    A_work( :, [ 2, 4 ] ) = A_work( :, [ 2, 4 ] ) * G3';
    A_work( [ 2, 3 ], : ) = G4 * A_work( [ 2, 3 ], : );
    A_work( :, [ 2, 3 ] ) = A_work( :, [ 2, 3 ] ) * G4';

    Q( :, [ 1, 2 ] ) = Q( :, [ 1, 2 ] ) * G1';
    Q( :, [ 1, 3 ] ) = Q( :, [ 1, 3 ] ) * G2';
    Q( :, [ 2, 4 ] ) = Q( :, [ 2, 4 ] ) * G3';
    Q( :, [ 2, 3 ] ) = Q( :, [ 2, 3 ] ) * G4';
else
    [ G1, H( [ 1, 2 ], 2 ) ] = planerot( H( [ 1, 2 ], 2 ) );
    H( [ 1, 2 ], 1 ) = G1 * H( [ 1, 2 ], 1 );
    [ G2, ~ ] = planerot( H( [ 1, 4 ], 2 ) );
    H( [ 1, 4 ], 1 ) = G2 * H( [ 1, 4 ], 1 );
    [ G3, H( [ 2, 3 ], 1 ) ] = planerot( H( [ 2, 3 ], 1 ) );
    [ G4, ~ ] = planerot( H( [ 2, 4 ], 1 ) );

    A_work( [ 1, 2 ], : ) = G1 * A_work( [ 1, 2 ], : );
    A_work( :, [ 1, 2 ] ) = A_work( :, [ 1, 2 ] ) * G1';
    A_work( [ 1, 4 ], : ) = G2 * A_work( [ 1, 4 ], : );
    A_work( :, [ 1, 4 ] ) = A_work( :, [ 1, 4 ] ) * G2';
    A_work( [ 2, 3 ], : ) = G3 * A_work( [ 2, 3 ], : );
    A_work( :, [ 2, 3 ] ) = A_work( :, [ 2, 3 ] ) * G3';
    A_work( [ 2, 4 ], : ) = G4 * A_work( [ 2, 4 ], : );
    A_work( :, [ 2, 4 ] ) = A_work( :, [ 2, 4 ] ) * G4';

    Q( :, [ 1, 2 ] ) = Q( :, [ 1, 2 ] ) * G1';
    Q( :, [ 1, 4 ] ) = Q( :, [ 1, 3 ] ) * G2';
    Q( :, [ 2, 3 ] ) = Q( :, [ 2, 3 ] ) * G3';
    Q( :, [ 2, 4 ] ) = Q( :, [ 2, 4 ] ) * G4';
end
% tolerance test
tol = max( max( abs( A_work( 3:4, 1:2 ) ) ) ) / norm( A );
thres = 10 * eps;
if  tol < thres
    A_work( 3:4, 1:2 ) = 0;
else 
    error( 'A fails to be reordered' );
end
end