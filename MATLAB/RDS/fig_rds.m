clear
log_start = -12; log_end = 12; num = 30;
k = 10.^(linspace(log_start, log_end, num)); 
g = 10.^(linspace(log_start, log_end, num));

random_num = 20;
lambda_r = zeros( num, 2 ); tol_r = zeros( num, 1 );
results_eig_r = zeros( num, num, 2 ); results_tol_r = zeros( num, num );

for i = 1:num 
    k_log = k( i ); 
    for j = 1:num 
        g_log = g( j );
        for m = 1:random_num
            [ A, eig_A ] = rand_test_matrices( k_log, g_log );
            [ A_work_r, ~, tol ] = swap22_sevp( A, 0 );            
            eig_r = eig_4_by_4( A_work_r ); 
                   
            lambda_r( m, 1 ) = abs( eig_A( 1 ) - eig_r( 2 ) ) / abs( eig_A( 1 ) );
            lambda_r( m, 2 ) = abs( eig_A( 2 ) - eig_r( 1 ) ) / abs( eig_A( 2 ) );
            tol_r( m ) = tol;
        end
        % Select the worst case 
        results_eig_r( num+1-i, j, : ) = max( lambda_r, [], 1 );
        results_tol_r( num+1-i, j ) = max( tol_r );

    end
end

%% Tol fig
data_tol = log10(results_tol_r);
log_min = floor( min( min(data_tol) ) ); log_max = ceil( max( max(data_tol) ) );

x_show_ticks = [.5, num/3+.5, 2*num/3+.5, num+.5];
y_show_ticks = x_show_ticks;
x_labels = {'$10^{-12}$', '$10^{-4}$', '$10^{4}$', '$10^{12}$'};
y_labels = {'$10^{12}$', '$10^{4}$', '$10^{-4}$', '$10^{-12}$'};

% Creat a new figure
figure('Position', [100, 100, 850, 650]);
imagesc(data_tol);
colormap(flipud(gray));

% Axis settings
ax = gca;
ax.XTick = x_show_ticks;      ax.YTick = y_show_ticks;
ax.XTickLabel = x_labels;     ax.YTickLabel = y_labels;
ax.TickLabelInterpreter = 'latex';
ax.TickDir = 'out';
ax.XGrid = 'off';             ax.YGrid = 'off';
ax.LineWidth = 2;
ax.FontSize = 16;

xlabel('$g$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$k$', 'Interpreter', 'latex', 'FontSize', 24);
title('RDS', 'FontSize', 20);

% Color range
clim([log_min, log_max]);

% Colorbar
c = colorbar;
tick_positions = unique([log_min:4:log_max, log_max]);
c.Ticks      = tick_positions;
c.TickLabels = arrayfun(@(v) sprintf('$10^{%d}$', v), ...
                            tick_positions, 'UniformOutput', false);
c.TickLabelInterpreter = 'latex';
c.FontSize   = 16;
c.LineWidth  = 2;
try c.TickDirection = 'out'; catch, end 

% Colorbar label
c.Label.String = sprintf( ...
        '$||\\check{A}_{2,1}||_{\\max} / ||A||_{\\max}$');
c.Label.Interpreter = 'latex';
c.Label.FontSize    = 20;
c.Label.Rotation    = 270;

axis image;
set(gcf, 'Color', 'white');

% Save current figure as PDF (filename includes ind to avoid overwriting)
exportgraphics(gcf, sprintf('RDS_tol.pdf'), 'ContentType', 'vector');

%% Eig fig
log_max_arr = [6, -2];
log_min = -16;

x_show_ticks = [.5, num/3+.5, 2*num/3+.5, num+.5];
y_show_ticks = x_show_ticks;
x_labels = {'$10^{-12}$', '$10^{-4}$', '$10^{4}$', '$10^{12}$'};
y_labels = {'$10^{12}$', '$10^{4}$', '$10^{-4}$', '$10^{-12}$'};

for ind = 1:2
    data_eig = log10( results_eig_r(:,:,ind) );

    % Create a new figure
    figure('Position', [100, 100, 850, 650]);
    imagesc(data_eig);
    colormap(flipud(gray));

    % Axis settings
    ax = gca;
    ax.XTick = x_show_ticks;      ax.YTick = y_show_ticks;
    ax.XTickLabel = x_labels;     ax.YTickLabel = y_labels;
    ax.TickLabelInterpreter = 'latex';
    ax.TickDir = 'out';
    ax.XGrid = 'off';             ax.YGrid = 'off';
    ax.LineWidth = 2;
    ax.FontSize = 16;

    xlabel('$g$', 'Interpreter', 'latex', 'FontSize', 24);
    ylabel('$k$', 'Interpreter', 'latex', 'FontSize', 24);
    title('RDS', 'FontSize', 20);

    % Color range
    log_max = log_max_arr(ind);
    clim([log_min, log_max]);

    % Colorbar
    c = colorbar;
    tick_positions = log_min:2:log_max;
    c.Ticks = tick_positions;
    c.TickLabels = arrayfun(@(v) sprintf('$10^{%d}$', v), ...
                            tick_positions, 'UniformOutput', false);
    c.TickLabelInterpreter = 'latex';
    c.FontSize = 16;
    c.LineWidth = 2;
    try c.TickDirection = 'out'; catch, end 

    % Colorbar label
    c.Label.String = sprintf( ...
        '$|\\lambda_{%d} - \\check{\\lambda}_{%d}| / |\\lambda_{%d}|$', ...
        ind, ind, ind);
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 20;
    c.Label.Rotation = 270;

    axis image;
    set(gcf, 'Color', 'white');

    % Save current figure as PDF (filename includes ind to avoid overwriting)
    exportgraphics(gcf, sprintf('MDS_eig_ind%d.pdf', ind), 'ContentType', 'vector');
end

%%
function  [A, eig_A] = rand_test_matrices(k, g) 
% Generate 4-by-4 block upper triangular having 2-by-2 blocks with complex-conjugate
% eigenvalues.
%
% Input 
% -----
%   k, g  :  parameters control the gap between the eigenvalues
% -----
% Output 
% -----
%   A  :  generated matrix
%   eig_A  :  eigenvalues of A
    a = randn(1); b = randn(1); r1 = randn(1); r2 = randn(1);
    A12 = randn(2);
    A11 = [a b*k; -b/k a];
    A22 = [a+r1*g (b+r2*g)*k; -(b+r2*g)/k a+r1*g];
    A = [A11 A12; zeros(2) A22];
    eig_A = eig_4_by_4(A); 
end

%%
function arr = eig_4_by_4( A )
% Compute and select eigenvalues of A with pisitive imaginary part.
% 
% Input
% -----
%   A  :  4-by-4 block upper triangular having 2-by-2 blocks with complex-conjugate
%         eigenvalues.
%
% Output
% ------
%   arr  :  2-by-1 vector with  .
if any( A( 3:4, 1:2 ) ~= 0, 'all')
    error( 'A( 3:4, 1:2 ) must be all zeros' );
end

eig1 = eig( A( 1:2, 1:2 ) );
eig2 = eig( A( 3:4, 3:4 ) );

if imag( eig1(1) ) >= 0
    val1 = eig1( 1 );
else
    val1 = eig1( 2 );  
end

if imag( eig2( 1 ) ) >= 0
    val2 = eig2( 1 );
else 
    val2 = eig2( 2 );  
end

arr = [ val1; val2 ];
end