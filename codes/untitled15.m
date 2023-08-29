A = imread("img.jpg");
Subplot (1,2,1);
imshow(A)
k= input ("enter k value" );
   if   A > k
        A =1;
  else 
       A= 0;
  end
subplot(1,2,2);
imshow(A);
