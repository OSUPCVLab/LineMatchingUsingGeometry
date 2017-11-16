clc
clear



im1 = imread('fig7.9c.gif');
im2 = imread('fig7.9d.gif');

[H,n] = homogtrans_ransac(im1,im2)