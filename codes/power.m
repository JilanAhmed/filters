figure;
c=1;
originalImage=imread('img.jpg');
gray_image = rgb2gray(originalImage);
double_value = im2double(gray_image);
out1= c*(double_value.^2); 
out2= c*(double_value.^0.4); 
out3= c*(double_value.^1.5); 

subplot(2,2,1);
imshow(gray_image);
title('original image');

subplot(2,2,2);
imshow((out1),[]);
title('(power=2)');

subplot(2,2,3);
imshow((out2),[]);
title('(power=0.4)');

subplot(2,2,4);
imshow((out3),[]);
title('(power=1.5)');
