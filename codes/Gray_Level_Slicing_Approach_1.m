%Gray Level Slicing Approach 1
figure
originalImage = imread("img.jpg");   
subplot(1, 2, 1),
  
% displaying the image
imshow(originalImage);
title("Original image");
gray_image = rgb2gray(originalImage);  
newImage = gray_image;
[rows cols] = size(gray_image);
for row_index=1:1:rows
    for col_index=1:1:cols
        if(gray_image(row_index,col_index)>=100 && gray_image(row_index,col_index)<=150)
            newImage(row_index,col_index) = 255;
        else
             newImage(row_index,col_index) = 0;
        end
    end
end
subplot(1, 2, 2),
imshow(newImage);
title("Gray Level Slicing Image")