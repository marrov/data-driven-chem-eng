%% Exercise 1:
% Understanding t-SNE

clc; clear; close all;

%% Read MNIST database

image_file = '../data/t10k-images.idx3-ubyte';
label_file = '../data/t10k-labels.idx1-ubyte';
[X, L] = read_mnist(image_file, label_file );

%% Reduce MNIST dataset size to reduce runtime

X = X(1:round(length(X)/3),:);
L = L(1:round(length(L)/3),:);

%% Run t-SNE

tic
rng default % for reproducibility
Y = tsne(X,'Algorithm','barneshut','NumPCAComponents',50);
toc

%% Plot t-SME with and without labels

figure(1)
gscatter(Y(:,1),Y(:,2))
figure(2)
numGroups = length(unique(L));
clr = hsv(numGroups);
gscatter(Y(:,1),Y(:,2),L,clr)