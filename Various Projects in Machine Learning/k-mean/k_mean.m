% clear all;
% close all;
matrix_cell_array = {train0, train1, train2, train3, train4, train5, train6, train7, train8, train9};
training_data = vertcat(matrix_cell_array{:});
training_data=double(training_data) / 255; 
%parameters
k = 10;
max_iter = 15; 
[num_points, num_features] = size(training_data);
centroids = rand(k, num_features);
labels = zeros(1, num_points); 
cost = zeros(1, max_iter); 

%k-means
for iter = 1:max_iter
    for i = 1:num_points
        min_dist = inf;
        for j = 1:k
            dist = sum((training_data(i,:) - centroids(j,:)).^2);
            if dist < min_dist
                min_dist = dist;
                labels(i) = j - 1;
            end
        end
    end
    
    cluster_sizes = histcounts(labels, 0:k);
    for j = 1:k
        cluster_points = training_data(labels == j-1, :);
        if ~isempty(cluster_points)
            centroids(j,:) = sum(cluster_points) / cluster_sizes(j);
        end
    end
    
    round_cost = 0;
    for j = 1:k
        cluster_points = training_data(labels == j-1, :);
        for i = 1:size(cluster_points, 1)
            diff = cluster_points(i,:) - centroids(j,:);
            round_cost = round_cost + sum(diff.^2);
        end
    end
    cost(iter) = round_cost;
end

figure;
plot(1:max_iter, cost, '-');
xlabel('Iterations', 'FontSize', 12);
ylabel('Cost', 'FontSize', 12);
title('Cost over number of iterations', 'FontSize', 15);

grid on;

