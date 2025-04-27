clear; close all; clc;

src = im2double(imread('../../srcImages/hazy3.png'));

A = airlightEstimation(src);
