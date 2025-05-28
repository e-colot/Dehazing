clear;
close all;

Hazy_true = imread("hazy1.png");
Hazy_image = im2double(Hazy_true);
Airlight_vector = AirlightDirection(Hazy_image);

Hazy_true = imread("hazy2.png");
Hazy_image = im2double(Hazy_true);
Airlight_vector_2 = AirlightDirection(Hazy_image);

Hazy_true = imread("hazy3.png");
Hazy_image = im2double(Hazy_true);
Airlight_vector_3 = AirlightDirection(Hazy_image);