figure
originalImage=imread('img.jpg');
gray_image = rgb2gray(originalImage)
imshow(img)

for i=1:row
        for j=1:col
           if  img(i,j)>128
                img(i,j)=255;
           else
                img(i,j)=0;
           end
        end
end 
figure
imshow(img);
