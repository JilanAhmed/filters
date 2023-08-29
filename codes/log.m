figure;
originalImage=imread("img.jpg");
gray_image = rgb2gray(originalImage);
double_value = im2double(gray_image);

out1= 2*log(1+double_value);
out2= 4*log(1+double_value);
out3= 6*log(1+double_value);

subplot(2,2,1);
imshow(gray_image);
title 'Original Image';

subplot(2,2,2);
imshow(out1);
title 'c=2';

subplot(2,2,3);
imshow(out2);
title 'c=4';

subplot(2,2,4);
imshow(out3);
title 'c=6';