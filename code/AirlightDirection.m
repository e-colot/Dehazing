function A_hat = AirlightDirection(I)
    % I: input RGB image (HxWx3, assumed double in [0, 1])
    disp('%%%%%%%%%%%%%%%%%BEGIN VECTOR DETERMINATION%%%%%%%%%%%%%%%%')
    % Parameters
    patch_size = 10;
    min_valid_patches = 30;
    %initial_patch_limit = 100;
    angle_threshold_deg = 15;

    % Step 1: Extract all non-edge patches
    gray = rgb2gray(I);
    edges = edge(gray, 'Canny');
    [H, W, ~] = size(I);
    valid_patches = struct('center', {}, 'direction', {}, 'eig1', {}, 'eig2', {}, 'dist_to_origin', {});
    patch_centers = [];

    % Precompute all 10x10 patches
    for i = 1:patch_size:H - patch_size + 1
        for j = 1:patch_size:W - patch_size + 1
            patch = I(i:i + patch_size - 1, j:j + patch_size - 1, :);
            patch_vec = reshape(patch, [], 3);

            % Skip patches with any edge pixels
            edge_patch = edges(i:i + patch_size - 1, j:j + patch_size - 1);
            if any(edge_patch(:))
                continue;
            end

            % Step 2: Robust PCA — first pass
            mu = mean(patch_vec, 1);
            centered = patch_vec - mu;
            [~, S, V] = svd(centered, 'econ');
            eigvals = diag(S).^2;

            % Step 3: Recompute PCA after removing top 20% farthest pixels
            dists = vecnorm(centered, 2, 2);
            keep = dists <= prctile(dists, 80);
            centered = centered(keep, :);
            if size(centered, 1) < 5
                continue;
            end
            [~, S2, V2] = svd(centered, 'econ');
            eigvals = diag(S2).^2;

            % PCA constraints
            principal = V2(:, 1);
            if any(principal < 0) && any(principal > 0)
                continue;
            end

            % Record patch info for thresholding
            valid_patches(end + 1).center = mu;
            valid_patches(end).direction = principal;
            valid_patches(end).eig1 = eigvals(1);
            valid_patches(end).eig2 = eigvals(2);
            valid_patches(end).dist_to_origin = norm(cross(mu, principal)) / norm(principal);
        end
    end

    fprintf('Number of valid patches found: %d\n', length(valid_patches));
    if isempty(valid_patches)
        warning('No valid patches found. Check edge detection or PCA constraints.');
        A_hat = [0;0;0];
        return;
    end

    % Threshold selection (initial filter)
    eig1_vals = [valid_patches.eig1];
    eig_ratio = [valid_patches.eig1] ./ max([valid_patches.eig2], 1e-5);
    dists = [valid_patches.dist_to_origin];

    % Normalize metrics
    eig1_norm = eig1_vals / max(eig1_vals);
    eig_ratio_norm = eig_ratio / max(eig_ratio);
    dists_norm = dists / max(dists);
    
    % Compute combined score (L2 norm)
    scores = sqrt(eig1_norm.^2 + eig_ratio_norm.^2 + dists_norm.^2);
    
    % Sort scores descending
    [sorted_scores, idx_scores] = sort(scores, 'descend');
    
    % Get the index of 50th best patch (or less if fewer patches)
    %desired_num = min(initial_patch_limit, length(scores));  % NOT WORKING
    % TO REACH 50 PATCHES MIN. -> DECIDED TO GO FOR ALL THE PATCHES
    desired_num = length(scores);
    threshold_idx = idx_scores(desired_num);
    
    % Set thresholds based on the patch at threshold_idx
    tau1 = eig1_vals(threshold_idx);
    tau2 = eig_ratio(threshold_idx);
    tau3 = dists(threshold_idx);

    % Refine patch list with combined thresholding
    filtered_patches = struct('center', {}, 'direction', {}, 'eig1', {}, 'eig2', {}, 'dist_to_origin', {});
    for i = 1:length(valid_patches)
        p = valid_patches(i);
        if p.eig1 >= tau1 && (p.eig1 / max(p.eig2, 1e-5)) >= tau2 && ...
           p.dist_to_origin >= tau3
            filtered_patches(end + 1) = p; %#ok<AGROW>
        end
    end

    fprintf('Number of filtered patches after thresholding: %d\n', length(filtered_patches));
    if isempty(filtered_patches)
        warning('No filtered patches remain after thresholding. Consider lowering thresholds.');
        A_hat = [0;0;0];
        return;
    end

    iter = 0;
    % Reduce until we get at least min_valid_patches
    while length(filtered_patches) > min_valid_patches && ...
        tau1 > 0 && tau2 > 0 && tau3 > 0
        tau1 = tau1 * 1.001;
        tau2 = tau2 * 1.001;
        tau3 = tau3 * 1.001;
        filtered_patches = struct('center', {}, 'direction', {}, 'eig1', {}, 'eig2', {}, 'dist_to_origin', {});
        for i = 1:length(valid_patches)
            p = valid_patches(i);
            if p.eig1 >= tau1 && (p.eig1 / max(p.eig2, 1e-5)) >= tau2 && ...
               p.dist_to_origin >= tau3
                filtered_patches(end + 1) = p; %#ok<AGROW>
            end
        end
    end

    fprintf('Number of filtered patches after reducing: %d\n', length(filtered_patches));
    if isempty(filtered_patches)
        warning('No filtered patches remain after thresholding. Consider lowering thresholds.');
        A_hat = [0;0;0];
        return;
    end

    A_hat_check = false;
    while A_hat_check == false
        % Remove redundant patches with similar directions
        final_patches = struct('center', {}, 'direction', {}, 'eig1', {}, 'eig2', {}, 'dist_to_origin', {});
        for i = 1:length(filtered_patches)
            if isempty(final_patches)
                final_patches = filtered_patches(i);
            else
                is_similar = false;
                for k = 1:length(final_patches)
                    angle = acosd(dot(filtered_patches(i).direction, final_patches(k).direction) / ...
                                  (norm(filtered_patches(i).direction) * norm(final_patches(k).direction)));
                    %fprintf('Angle between patch %d and %d: %.2f degrees\n', i, k, angle);
                    if angle <= angle_threshold_deg
                        is_similar = true;
                        break;
                    end
                end
    
                if ~is_similar
                    final_patches(end + 1) = filtered_patches(i); %#ok<AGROW>
                end
            end
        end
        fprintf('Number of final patches after direction redundancy removal: %d\n', length(final_patches));
        if length(final_patches) < 2
            warning('Not enough final patches to estimate airlight direction.');
            A_hat = [0;0;0];
            return;
        end
    
        % Estimate candidate A directions via plane intersections
        A_candidates = [];
        for i = 1:length(final_patches)
            for j = i+1:length(final_patches)
                n1 = cross(final_patches(i).center, final_patches(i).direction);
                n2 = cross(final_patches(j).center, final_patches(j).direction);
                A_dir = cross(n1, n2);
                if norm(A_dir) > 1e-6
                    A_candidates(:, end + 1) = A_dir / norm(A_dir); %#ok<AGROW>
                end
            end
        end
        fprintf('Number of candidate A directions computed: %d\n', size(A_candidates, 2));
        if isempty(A_candidates)
            warning('No candidate airlight directions found from patch intersections.');
            A_hat = [0;0;0];
            return;
        end
    
        % Choose best A_dir by lowest median distance to patch lines
        min_median = Inf;
        A_hat = [0; 0; 0];
        for j = 1:size(A_candidates, 2)
            candidate = A_candidates(:, j);
            dists = zeros(length(final_patches), 1);
            for i = 1:length(final_patches)
                % Compute shortest distance between candidate ray and patch line
                u = candidate;
                v = final_patches(i).direction;
                p1 = zeros(3, 1);
                p2 = final_patches(i).center';
                cross_uv = cross(u, v);
                w = p2 - p1;
                proj = dot(w, cross_uv) * cross_uv / (norm(cross_uv)^2);
                d = norm(w - proj);
                dists(i) = d;
            end
            med = median(dists);
            if med < min_median
                min_median = med;
                A_hat = candidate;
            end
        end
        if any(A_hat < 0) && any(A_hat > 0)
            angle_threshold_deg = angle_threshold_deg + 5;
        elseif size(A_candidates, 2) > 10
            angle_threshold_deg = angle_threshold_deg + 5;
        else
            break;
        end
    end
    
    % Case if all A_hat components are negatives
    if any(A_hat < 0)
        A_hat = -A_hat;
    end

    disp('%%%%%%%%%%%%%%%%%END VECTOR DETERMINATION%%%%%%%%%%%%%%%%')

    figure;
    hold on;
    grid on;
    xlabel('R');
    ylabel('G');
    zlabel('B');
    title('Patch Lines in RGB Space and Estimated Airlight Direction');
    
    % Plot all selected patch lines in blue (thicker lines)
    for i = 1:length(final_patches)
        p = final_patches(i);
        line_pts = [p.center' - 0.2 * p.direction, p.center' + 0.2 * p.direction];
        plot3(line_pts(1,:), line_pts(2,:), line_pts(3,:), 'b-', 'LineWidth', 2);
    end
    
    % Plot the estimated A_hat direction (ray from origin) in red
    quiver3(0, 0, 0, A_hat(1), A_hat(2), A_hat(3), ...
        0.5, 'r', 'LineWidth', 3, 'MaxHeadSize', 1);
    
    legend('Patch RGB lines', 'Estimated Â');
    axis equal;
    xlim([0 1]);
    ylim([0 1]);
    zlim([0 1]);
    view(3); % ensures 3D view
end

