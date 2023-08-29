
figure
originalImage=imread('img.jpg');
a=double(a);
s=size(a);
subplot(1,2,1)
imshow(a/255);

for n=1:s(1,1)
    for m=1:s(1,2)
        ot(n,m)=t(a(n,m)+1);
    end
end

subplot(1,2,2)
imshow(ot./255)


