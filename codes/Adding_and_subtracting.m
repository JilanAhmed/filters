% Adding and subtracting
figure;
originalImage = imread("img.jpg");   
subplot(3,1, 1),
% displaying the image
imshow(originalImage);
title("Original image");
gray_image = rgb2gray(originalImage); 
AddImage = gray_image+120;
SubtractImage = gray_image-120;
subplot(3,1,2),
imshow(AddImage);
title("Image after the add process")
subplot(3,1,3),
imshow(SubtractImage);
title("Image after the subtract process")