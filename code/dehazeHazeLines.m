function [J, t] = dehazeHazeLines(I, A)
% Dehazing using Haze-Lines (Berman et al., CVPR 2016)
% I: input hazy image (HxWx3)
% A: airlight (1x3)
% J: dehazed image
% t: estimated transmission map

[h,w,~] = size(I);

% Step 1: Center on A
    IA = I - reshape(A,1,1,3);
    IA = reshape(IA, h*w, 3); % Reshape to 2D array for processing

% Step 2: Calculate radius
    R = max(sqrt(sum(IA.^2,2)), 1e-6); % No division by zero

% Step 3: Define haze lines by clustering the pixels
    unitaryRadius = reshape(IA./ R, [], 3); % Normalize

    n_points = 1000;

    % load pre-calculated uniform tesselation of the unit-sphere
    fid = fopen(['TR',num2str(n_points),'.txt']);
    points = cell2mat(textscan(fid,'%f %f %f')) ;
    fclose(fid);

    % build a KD-tree (fast for nearest neighbor search)
    mdl = KDTreeSearcher(points);
    % find the nearest neighbors in the KD-tree
    idx = knnsearch(mdl, unitaryRadius);

% Step 4-7: Estimate Initial Transmission
    R_max = accumarray(idx, R(:), [n_points, 1], @max); % maximal R in each haze-line -> length(R_max) = 1000
    t_hat = R ./ R_max(idx);
    % limit t_hat to [0.1, 1]
    t_hat = min(max(t_hat, 0.1), 1);
    test = t_hat;

% Step 8: Regularization
    % lower bound correction
    t_LB = reshape(I, h*w, 3) ./ reshape(A, 1, 3);
    t_LB = 1 - min(t_LB, [], 2); % lower bound for each pixel

    t_LB_hat = max(t_hat, t_LB);

    % code taken as such from https://github.com/danaberman/non-local-dehazing
    % comments added for clarity
    bin_count       = accumarray(idx, 1, [n_points, 1]);       % number of pixels per haze-line
    bin_count_map   = reshape(bin_count(idx), h, w);
    bin_eval_fun    = @(x) min(1, x/50);                    % set a threshold for small haze-lines
    % if the haze-line has less than 50 pixels, it is not reliable and should be disregarded
                        
    % Calculate standard deviation of radius in each haze-line
    K_std = accumarray(idx, R(:), [n_points, 1], @std);
    radius_std = reshape(K_std(idx), h, w);

    % remove pixels that have a small standard deviation in radius
    radius_eval_fun = @(r) min(1, 3*max(0.001, r-0.1));
    radius_reliability = radius_eval_fun(radius_std./max(radius_std(:)));

    % combine conditions of reliability based on bin count and radius standard deviation
    data_term_weight   = bin_eval_fun(bin_count_map).*radius_reliability;
    lambda = 0.1;

    % reshape t_LB_hat to match the size of I
    t_LB_hat = reshape(t_LB_hat, h, w);

    % solve optimization problem (Eq. (15))
    t = wls_optimization(t_LB_hat, data_term_weight, I, lambda);

    % limit t to [0.1, 1]
    t = min(max(t, 0.1), 1);

    t = reshape(test, h, w);

J = (I - (1 - t) .* reshape(A, 1, 1, 3)) ./ t;

end