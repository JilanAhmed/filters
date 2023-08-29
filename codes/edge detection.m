
%% Laplacian filter
 originalImage = imread("img.jpg");  
 % gray_image = rgb2gray(originalImage);  
 subplot(1, 2, 1),
 % displaying the image
 imshow(originalImage);
 title("Original image");
 
 gray_image = double(originalImage);
 [rows,cols]=size(gray_image);
 mask = [0,-1,0;-1,5,-1;0,-1,0];
 out = gray_image;
 for i=2:rows-1
  for j=2:cols-1
      temp = mask.*gray_image(i-1:i+1,j-1:j+1);
      value = sum(temp(:));
      out(i, j)= value;
 end
 end
 out = uint8(out);
 subplot(1, 2, 2),
 imshow(out);
 title("sharp image");
%% Sobel operator

%originalImage = imread("img.jpg");  
% gray_image = rgb2gray(originalImage);  
%subplot(1, 2, 1),
% displaying the image
%imshow(originalImage);
%title("Original image");

%gray_image = double(originalImage);
%[rows,cols]=size(gray_image);
% mask = [-1 -2 -1;0 0 0;1 2 1];
%mask = [-1 0 1;-2 0 2;-1 0 1];
%out = gray_image;
%for i=2:rows-1
 %for j=2:cols-1
     %temp = mask.*gray_image(i-1:i+1,j-1:j+1);
     %value = sum(temp( ));
     %out(i, j)= value;
%end
%end
%out = uint8(out);
%subplot(1, 2, 2),
%imshow(out);
%title("edge detection");


