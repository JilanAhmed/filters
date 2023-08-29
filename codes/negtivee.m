figure;
originalImage = imread("img.jpg");   
subplot(1, 2, 1),
  
imshow(originalImage);
title("Original image");
  
L = 2 ^ 8;    
                  
negtiveImage = (L - 1) - originalImage;
subplot(1, 2, 2),
  
imshow(negtiveImage);
title("Negative Image")