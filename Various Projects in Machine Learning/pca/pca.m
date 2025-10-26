% clear all
% close all
% Divide the dataset into training and testing sets with a 9:2 ratio:
% Training uses the first 9 images out of 11 for each person, while testing uses the last 2.
num_people = 15;
num_train = 9;
num_test = 2;
data_dim = 1024;

train_data = zeros(num_people * num_train, data_dim);
labelTrain_data = zeros(num_people * num_train, 1);
test_data = zeros(num_people * num_test, data_dim);
labelTest_data = zeros(num_people * num_test, 1);

for person = 1:num_people
    start_idx = 1 + (person - 1) * 11;
    end_train_idx = start_idx + num_train - 1;
    end_test_idx = start_idx + 10;

    train_data(1 + (person - 1) * num_train : person * num_train, :) = faces(start_idx:end_train_idx, :);
    labelTrain_data(1 + (person - 1) * num_train : person * num_train) = labeles(start_idx:end_train_idx);

    test_data(1 + (person - 1) * num_test : person * num_test, :) = faces(end_train_idx + 1:end_test_idx, :);
    labelTest_data(1 + (person - 1) * num_test : person * num_test) = labeles(end_train_idx + 1:end_test_idx);
end

% Normalize the data to center around zero:
meanTrain = mean(train_data);
shiftedTrain = train_data - meanTrain;
meanTest = mean(test_data);
shiftedTest = test_data - meanTest;

% Perform PCA:
% finding eigenvalues and cov-matrixes:
covTrain = (1 / size(shiftedTrain, 1)) * (shiftedTrain' * shiftedTrain);
[V, ~] = eig(covTrain); 
V = fliplr(V); 

dim = 300; 
SuccessRate = zeros(1, dim);

% PCA and nearest neighbor:
for d = 1:dim
    projTrain = shiftedTrain * V(:, 1:d); 
    projTest = shiftedTest * V(:, 1:d); 
    correct_guesses = 0;

    for test_idx = 1:size(projTest, 1)
        distances = sum((projTrain - projTest(test_idx, :)).^2, 2); 
        %NNE:
        [~, min_idx] = min(distances);
        %finding success rate comparing labels:
        if labelTrain_data(min_idx) == labelTest_data(test_idx)
            correct_guesses = correct_guesses + 1;
        end
    end

    SuccessRate(d) = correct_guesses / size(projTest, 1); 
end

% Visualize original dataset images:
Pic0 = zeros(32 * num_people, 32 * 11);
for person = 1:num_people
    for img_idx = 1:11
        img = reshape(faces(img_idx + (person - 1) * 11, :), 32, 32);
        Pic0(1 + 32 * (person - 1):32 * person, 1 + 32 * (img_idx - 1):32 * img_idx) = img;
    end
end
figure(1);
imshow(Pic0, []);

% Visualize mean-centered training images:
Pic1 = zeros(32 * num_people, 32 * num_train);
for person = 1:num_people
    for img_idx = 1:num_train
        img = reshape(shiftedTrain(img_idx + (person - 1) * num_train, :), 32, 32);
        Pic1(1 + 32 * (person - 1):32 * person, 1 + 32 * (img_idx - 1):32 * img_idx) = img;
    end
end
figure(2);
imshow(Pic1, []);

% Plot success rate versus number of dimensions:
figure(3);
plot(1:dim, SuccessRate, 'LineWidth', 1.5);
title('Success rate \propto Number of dimensions');
xlabel('Number of dimensions');
ylabel('Success rate');
legend('Stable with ~120 dimensions', 'Location', 'east');

