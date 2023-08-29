%% Bit Plan Slicing
figure;
originalImage = imread("img.jpg");   
subplot(2, 5, 1),
% displaying the image
imshow(originalImage);
title("Original image");
gray_image = rgb2gray(originalImage); 
[rows cols] = size(gray_image);
newImage = zeros(rows,cols,8);
for k=1:8
    for row_index=1:1:rows
        for col_index=1:1:cols
            newImage(row_index,col_index,k)=bitget(gray_image(row_index,col_index),k);
        end
    end
subplot(2, 5, k+1),
imshow(newImage(:,:,k));
title("Image for bit number "+k)
end