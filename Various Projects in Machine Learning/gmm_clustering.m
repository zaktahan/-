function gmm_clustering()
    rng(42); 
    %parameters
    alpha = [0.45, 0.30, 0.25]; 
    mu = [2, 3; 6, 3; 3, 5]; 
    sigma(:,:,1) = [0.6, 0; 0, 0.4]; 
    sigma(:,:,2) = [0.8, 0; 0, 0.6]; 
    sigma(:,:,3) = [0.7, 0; 0, 0.5]; 
    [data, true_labels] = generate_data(alpha, mu, sigma, 1000);
    manual_kmeans(data, 3, 100, true_labels);
    manual_gmm(data, 3, 100, true_labels);
end
function [data, labels] = generate_data(alpha, mu, sigma, n)
    k = length(alpha); 
    labels = zeros(n, 1);
    data = zeros(n, size(mu, 2));
    for i = 1:n
        z = find(rand <= cumsum(alpha), 1); 
        labels(i) = z;
        data(i, :) = mvnrnd(mu(z, :), sigma(:, :, z), 1);
    end
end
%k means algorithm
function manual_kmeans(data, k, max_iter, true_labels)
    [n, d] = size(data);
    centroids = data(randperm(n, k), :);
    
    labels = zeros(n, 1);
    
    for iter = 1:max_iter
        for i = 1:n
            distances = sum((data(i, :) - centroids).^2, 2);
            [~, labels(i)] = min(distances);
        end
        for j = 1:k
            cluster_points = data(labels == j, :);
            if ~isempty(cluster_points)
                centroids(j, :) = mean(cluster_points, 1);
            end
        end
    end
    accuracy = calculate_accuracy(labels, true_labels);
    fprintf('K-means Accuracy: %.2f%%\n', accuracy * 100);
    
    plot_clusters(data, labels, 'Manual K-means Clustering');
end
%mixture of gaussian algorithm
function manual_gmm(data, k, max_iter, true_labels)
    [n, d] = size(data);
    
    phi = ones(1, k) / k; 
    mu = data(randperm(n, k), :); 
    sigma = repmat(eye(d), [1, 1, k]); 
    responsibilities = zeros(n, k);
    
    for iter = 1:max_iter
        for j = 1:k
            responsibilities(:, j) = phi(j) * mvnpdf(data, mu(j, :), sigma(:, :, j));
        end
        responsibilities = responsibilities ./ sum(responsibilities, 2);
        
        for j = 1:k
            wj = sum(responsibilities(:, j));
            phi(j) = wj / n;
            mu(j, :) = sum(responsibilities(:, j) .* data) / wj;
            
            sigma(:, :, j) = zeros(d, d);
            for i = 1:n
                diff = data(i, :) - mu(j, :);
                sigma(:, :, j) = sigma(:, :, j) + responsibilities(i, j) * (diff' * diff);
            end
            sigma(:, :, j) = sigma(:, :, j) / wj;
        end
    end
    
    [~, labels] = max(responsibilities, [], 2);
    
    accuracy = calculate_accuracy(labels, true_labels);
    fprintf('GMM Accuracy: %.2f%%\n', accuracy * 100);
    
    plot_clusters(data, labels, 'Manual GMM Clustering');
end
%finding accuracy
function accuracy = calculate_accuracy(predicted_labels, true_labels)
    k = max(true_labels);
    mapping = zeros(1, k);
    for i = 1:k
        most_common = mode(true_labels(predicted_labels == i));
        mapping(i) = most_common;
    end
    
    adjusted_labels = arrayfun(@(x) mapping(x), predicted_labels);
    
    accuracy = sum(adjusted_labels == true_labels) / length(true_labels);
end

function plot_clusters(data, labels, title_str)
    figure;
    gscatter(data(:, 1), data(:, 2), labels);
    title(title_str);
    xlabel('x1'); ylabel('x2');
    grid on;
end