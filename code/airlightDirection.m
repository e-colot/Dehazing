function A = airlightDirection(img)

    dim = size(img);

    centroids = zeros(0, 3); % centroids of the patches
    PC = zeros(0, 3); % principal components of the patches
    % treshold values for each saved patch
    crit1 = [];
    crit2 = [];
    crit3 = [];

    edges = edge(mean(img, 3), "canny");
    patchSize = 10;

    for i = 1:floor(dim(1)/patchSize)
        for j = 1:floor(dim(2)/patchSize)

            patch = img(1+patchSize*(i-1):patchSize*i,1+patchSize*(j-1):patchSize*j,:);
            centroid = mean(mean(patch, 1), 2);
            centroid = reshape(centroid, 1, 3);

            % PCA (Principal Component Analysis)
            patch_reshaped = reshape(patch, [], 3);
            cov_matrix = cov(double(patch_reshaped));
            [eigenvectors, eigenvalues] = eig(cov_matrix);
            % Find the principal component (largest eigenvector)
            [~, max_idx] = max(diag(eigenvalues));
            principal_component = eigenvectors(:, max_idx);

            % Project pixels onto the principal component
            projections = patch_reshaped * principal_component;
            distances = abs(projections - mean(projections));

            % Sort distances and discard the farthest 20% pixels
            [~, sorted_indices] = sort(distances, 'ascend');
            cutoff_index = round(0.8 * length(sorted_indices));
            filtered_patch = patch_reshaped(sorted_indices(1:cutoff_index), :);

            % Recompute PCA with the filtered patch
            cov_matrix_filtered = cov(double(filtered_patch));
            [eigenvectors_filtered, eigenvalues_filtered] = eig(cov_matrix_filtered);
            [eigenvalues_filtered_sorted, sort_indices] = sort(diag(eigenvalues_filtered), 'descend');
            eigenvectors_filtered = eigenvectors_filtered(:, sort_indices);

            principal_component = eigenvectors_filtered(:, 1);
            if sum(principal_component >= 0) ~= 3
                % 'Positive principal component' rule broken
                continue;
            end

            numberOfEdges = sum(edges(1+patchSize*(i-1):patchSize*i,1+patchSize*(j-1):patchSize*j), 'all');
            if numberOfEdges > 0
                % 'Patches do not contain an edge' rule broken
                continue;
            end

            centroids = [centroids; centroid];
            PC = [PC; principal_component'];

            crit1 = [crit1; eigenvalues_filtered_sorted(1)];
            crit2 = [crit2; eigenvalues_filtered_sorted(1)/eigenvalues_filtered_sorted(2)];
            d = norm(cross(centroid, principal_component)) / norm(principal_component);
            crit3 = [crit3; d];

        end
    end

    crit1Sorted = sort(crit1, 'descend');
    crit2Sorted = sort(crit2, 'descend');
    crit3Sorted = sort(crit3, 'descend');

    thresholds = [crit1Sorted(51), crit2Sorted(51), crit3Sorted(51)];
    validIndexes = [];
    
    for i = 1:size(centroids, 1)
        if crit1(i) < thresholds(1) || crit2(i) < thresholds(2) || crit3(i) < thresholds(3)
            continue;
        end
        validIndexes = [validIndexes; i];
    end

    while size(validIndexes, 1) < 10
        thresholds = thresholds * 0.97;
        validIndexes = [];
        for i = 1:size(centroids, 1)
            if crit1(i) < thresholds(1) || crit2(i) < thresholds(2) || crit3(i) < thresholds(3)
                continue;
            end
            validIndexes = [validIndexes; i];
        end
    end

    % an angle criterion is described in the reference paper
    % as the line is 3 dimensional, the condition is quite complex to implement

    figure;
    hold on;
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('3D Lines through Centroids with Principal Component Directions');

    for i = 1:length(validIndexes)
        idx = validIndexes(i);
        point = centroids(idx, :);
        direction = PC(idx, :);

        % Define the line in both directions
        t = -10:0.1:10; % Parameter for the line
        line_points = point + t' * direction;

        % Plot the line
        plot3(line_points(:, 1), line_points(:, 2), line_points(:, 3));
    end

    xlim([0 1]);
    ylim([0 1]);
    zlim([0 1]);
    hold off;

    % Compute candidate orientation vectors by intersecting planes
    candidateOrientations = [];
    for i = 1:length(validIndexes)
        for j = i+1:length(validIndexes)
            idx1 = validIndexes(i);
            idx2 = validIndexes(j);

            % Centroids and principal components of the two patches
            point1 = centroids(idx1, :);
            direction1 = PC(idx1, :);
            point2 = centroids(idx2, :);
            direction2 = PC(idx2, :);

            % Compute the normal vector of the plane defined by each line
            normal1 = cross(direction1, point1);
            normal2 = cross(direction2, point2);

            % Compute the candidate orientation vector as the intersection of the planes
            candidateOrientation = cross(normal1, normal2);
            if norm(candidateOrientation) > 0
                candidateOrientation = candidateOrientation / norm(candidateOrientation);
                candidateOrientations = [candidateOrientations; candidateOrientation];
            end
        end
    end

    % Compute distances between patch lines and candidate rays
    distances = zeros(size(centroids, 1), size(candidateOrientations, 1));
    for i = 1:size(centroids, 1)
        point = centroids(i, :);
        direction = PC(i, :);

        for j = 1:size(candidateOrientations, 1)
            rayDirection = candidateOrientations(j, :);

            % Compute the nearest points on the line and the ray
            t = dot(rayDirection, (point - dot(point, rayDirection) * rayDirection)) / dot(rayDirection, rayDirection);
            nearestPointOnRay = t * rayDirection;
            nearestPointOnLine = point + dot(nearestPointOnRay - point, direction) * direction;

            % Compute the Euclidean distance
            distances(i, j) = norm(nearestPointOnLine - nearestPointOnRay);
        end
    end

    % Select the ray with the lowest median distance
    medianDistances = median(distances, 1);
    [~, bestRayIndex] = min(medianDistances);
    A = candidateOrientations(bestRayIndex, :);

    hold on;
    quiver3(0, 0, 0, A(1), A(2), A(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    hold off;

end
