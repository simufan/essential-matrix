figure(1);
im1 = imread('0000.scaled.jpg');
imshow(im1);
hold on;
plot(m1(:,1), m1(:,2), 'R+', 'LineWidth', 2, 'MarkerSize',10);
hold off;


figure(2);
im2 = imread('0001.scaled.jpg');
imshow(im2);
hold on;
plot(m2(:,1), m2(:,2), 'R+', 'LineWidth', 2, 'MarkerSize',10);
hold off;

k = [2759.48 0 1520.69; 
0 2764.16 1006.81;
0 0 1];

Width = 3072; %image width
Height = 2048; %image height

