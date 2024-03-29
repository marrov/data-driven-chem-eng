%% Exercise 1:
% Understanding t-SNE

clc; clear; close all;

%% Read MNIST database

image_file = '../data/images.idx3-ubyte';
label_file = '../data/labels.idx1-ubyte';
[X, L] = read_mnist(image_file, label_file );

%% Reduce MNIST dataset size to reduce runtime

X = X(1:round(length(X)/2),:);
L = L(1:round(length(L)/2),:);

%% Run t-SNE

rng default % for reproducibility
Y = tsne(X,'Algorithm','barneshut','NumPCAComponents',50);

%% Plot t-SME with and without labels

figure(1)
gscatter(Y(:,1),Y(:,2))
figure(2)
numGroups = length(unique(L));
clr = hsv(numGroups);
gscatter(Y(:,1),Y(:,2),L,clr)

%% Try for youself...

% Different values of perplexity
% Y = tsne(X, ...

% Different distance metrics
% Y = tsne(X, ...

% At least one of the following: PCA components, learning rate, exageration
% Y = tsne(X, ...

% Not using "rng default"
% Y = tsne(X, ...