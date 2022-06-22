%% Exercise 1:
% Understanding t-SNE

clc; clear; close all;

image_file = 'data/t10k-images.idx3-ubyte';
label_file = 'data/t10k-labels.idx1-ubyte';

[X, L] = read_mnist(image_file, label_file );

rng default % for reproducibility
Y = tsne(X,'Algorithm','barneshut','NumPCAComponents',50);

figure
numGroups = length(unique(L));
clr = hsv(numGroups);
gscatter(Y(:,1),Y(:,2),L,clr)