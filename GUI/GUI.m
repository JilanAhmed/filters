function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 02-Mar-2022 19:33:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;
A = imread("img.jpg");
subplot (1,2,1);
imshow(A)
sz = size(A);
xg = 1:sz(1);
yg = 1:sz(2);
F = griddedInterpolant({xg,yg},double(A));
xq = (0:5/6:sz(1))';
yq = (0:5/6:sz(2))';
vq = uint8(F({xq,yq}));
imshow(vq);
title("Higher Resolution");
xq = (0:1.55:sz(1))';
yq = (0:1.55:sz(2))';
vq = uint8(F({xq,yq}));
subplot(1,2,2);
imshow(vq);
title("Lower Resolution")

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

img = imread("img.jpg");
k = 8;
figure
while (k > 0)
 target_levels = 2^k;
 target_compr_factor = 256 / target_levels;
 reduced_image = uint8(floor(double(img)/256 * target_levels) * target_compr_factor);
 subplot(3, 3, k); 
 imshow(reduced_image, [0 255]);
 if (k==1)
      title('Black & White');
 else
      title(['Grey-level resolution 2^',num2str(k)]);
 end
 k = k - 1;
end

% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
originalImage = imread("img.jpg");
[rows cols matricesNo] = size(originalImage);
SamplingFactor = 16;
for metricesIndex=1:1:matricesNo
    resizedImage(:,:,metricesIndex) = subSampling(originalImage(:,:,metricesIndex),SamplingFactor);
end
imshow(resizedImage);
imwrite(resizedImage,'resizedImage.png');
function outImage = subSampling(image, subSamplingFactor)
[rows cols] = size(image);
outImage = image(1:subSamplingFactor:rows,1:subSamplingFactor:cols);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
Img1 = imread('histo.png');
Img = rgb2gray(Img1);
subplot(2,1,1);
imshow(Img)
title("Original Image");

[row, col] = size(Img);
 for x=1:row 
        for y=1:col
           if  Img(x,y)>128
                Img(x,y)=255;
           else
                Img(x,y)=0;
           end
        end
end 
subplot(2,1,2);
imshow(Img)
title("Thresholing");


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
originalImage = imread("thre.png"); 
 gray_image = rgb2gray(originalImage);  
 subplot(1, 2, 1),

 imshow(gray_image);
 title("Original image");
 
 [rows,cols]=size(gray_image);
 counts=zeros(1,256);
 for i=1:rows
  for j=1:cols
     grayLevel=gray_image(i,j);
     counts(grayLevel+1)=counts(grayLevel+1)+1;
 end
 end

 subplot(1, 2, 2),
 grayLevels = 0 : 255;
 bar(grayLevels, counts, 'BarWidth', 1, 'FaceColor', 'b');
 xlabel('Gray Level', 'FontSize', 20);
 ylabel('Pixel Count', 'FontSize', 20);
 title('Histogram', 'FontSize', 20);
 grid on;





% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
originalImage = imread("sobel.png");  
 GIm = rgb2gray(originalImage);
 numofpixels=size(GIm,1)*size(GIm,2);
 subplot(1, 2, 1),
 imshow(GIm);
 title("Original image");
 
 HIm=uint8(zeros(size(GIm,1),size(GIm,2)));
 
 freq=zeros(256,1);
 
 probc=zeros(256,1);
 
 output=zeros(256,1);
 
 for i=1:size(GIm,1)

    for j=1:size(GIm,2)

        value=GIm(i,j);

        freq(value+1)=freq(value+1)+1;

     end
 
 end
 
 
 sum=0;
 no_bins=255;

 
 for i=1:size(freq)
 
    sum=sum+freq(i);
 
   probc(i)=sum/numofpixels;

   output(i)=round(probc(i)*no_bins);

 end
 
 for i=1:size(GIm,1)
 
     for j=1:size(GIm,2)

  HIm(i,j)=output(GIm(i,j)+1);

   end

 end
 subplot(1, 2, 2),
 imshow(HIm);
title('Histogram equalization');


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
gray_image = imread("img.jpg");    
subplot(1, 2, 1),

imshow(gray_image);
title("Original image");

[rows,cols]=size(gray_image);
out=gray_image;

for i=2:rows-1
 for j=2:cols-1
     temp = [gray_image(i-1, j-1) gray_image(i-1, j) gray_image(i-1, j + 1) gray_image(i, j-1) gray_image(i, j) gray_image(i, j + 1) gray_image(i + 1, j-1) gray_image(i + 1, j) gray_image(i + 1, j + 1)];
     temp = mean(temp);
     out(i, j)= temp;
 
end
end

subplot(1, 2, 2),
imshow(out);
title("Average Image");


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Gray Level Slicing Approach 2
figure;
originalImage = imread("log.png");   
subplot(1, 2, 1),
%displaying the image
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
             newImage(row_index,col_index) = gray_image(row_index,col_index);
        end
    end
end

subplot(1, 2, 2),
imshow(newImage);
title("Gray Level Slicing Image approach2")

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;
originalImage=imread("power.png");
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


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;
originalImage = imread("grayscale.png");   
subplot(1, 2, 1),
  
imshow(originalImage);
title("Original image");
  
L = 2 ^ 8;    
                  
negtiveImage = (L - 1) - originalImage;
subplot(1, 2, 2),
  
imshow(negtiveImage);
title("Negative Image")


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;
c=1;
o=imread('img.jpg');
double_value = im2double(o);
out1= c*(double_value.^2); 
out2= c*(double_value.^0.4); 
out3= c*(double_value.^1.5); 

subplot(2,2,1);
imshow(o);
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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure
originalImage = imread("laplacian.png");  
% gray_image = rgb2gray(originalImage);  
subplot(1, 2, 1),
% displaying the image
imshow(originalImage);
title("Original image");

gray_image = double(originalImage);
[rows,cols]=size(gray_image);
mask = [0,1,0;1,-4,1;0,1,0];
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


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Adding and subtracting
figure;
originalImage = imread("sharpen.png");   
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

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Bit Plan Slicing
figure;
originalImage = imread("neg.png");   
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

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
xrayImg = imread('median.jpg');
subplot(2,1,1);
imshow(xrayImg)
title("Original Image");

xrayImg=rgb2gray(xrayImg);

[row, col] = size(xrayImg);
rmin=min(min(xrayImg));
rmax=max(max(xrayImg));

for x=1:row 
        for y=1:col
           xrayImg(x,y)=((255)/(rmax-rmin))*((xrayImg(x,y))-(rmin));
        end
end
subplot(2,1,2);
imshow(xrayImg)
title("Contrast Image")



% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% ILPF
figure
input_image = imread("ih.png");  
[M, N] = size(input_image , [1 2]);
FT_img = fft2(im2double(input_image));
D0 = 15; 

% Designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;

% MATLAB library function meshgrid(v, u) returns
% 2D grid which contains the coordinates of vectors
% v and u. Matrix V with each row is a copy 
% of v, and matrix U with each column is a copy of u
[V, U] = meshgrid(v, u);
  
% Calculating Euclidean Distance
D = sqrt(U.^2+V.^2);
  
% Comparing with the cut-off frequency and 
% determining the filtering mask
H = double(D <= D0);
  
% Convolution between the Fourier Transformed
% image and the mask
G = H.*FT_img;
  
% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)  
output_image = real(ifft2(double(G)));
% Displaying Input Image and Output Image
subplot(2, 1, 1), imshow(input_image),
subplot(2, 1, 2), imshow(output_image, [ ]);

% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% BLPF
figure;
input_image = imread("bh.png");  
[M, N] = size(input_image , [1 2]);
FT_img = fft2(im2double(input_image));
D0 = 15; 
n=2*2;
% Designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;

% MATLAB library function meshgrid(v, u) returns
% 2D grid which contains the coordinates of vectors
% v and u. Matrix V with each row is a copy 
% of v, and matrix U with each column is a copy of u
[V, U] = meshgrid(v, u);
  
% Calculating Euclidean Distance
D = sqrt(U.^2+V.^2);

D = D./ D0;
% Comparing with the cut-off frequency and 
% determining the filtering mask
H = 1./((1+D).^n);
  
% Convolution between the Fourier Transformed
% image and the mask
G = H.*FT_img;
  
% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)  
output_image = real(ifft2(double(G)));
% Displaying Input Image and Output Image
subplot(2, 1, 1), imshow(input_image),
subplot(2, 1, 2), imshow(output_image, [ ]);

% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% GLPF
figure;
input_image = imread("gh.png");  
[M, N] = size(input_image , [1 2]);
FT_img = fft2(im2double(input_image));
D0 = 15; 
D0 = (D0^2)*2;
% Designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;

% MATLAB library function meshgrid(v, u) returns
% 2D grid which contains the coordinates of vectors
% v and u. Matrix V with each row is a copy 
% of v, and matrix U with each column is a copy of u
[V, U] = meshgrid(v, u);
  
% Calculating Euclidean Distance
D = sqrt(U.^2+V.^2);

D = -D.^2;
% Comparing with the cut-off frequency and 
% determining the filtering mask
H = exp(D/D0);
  
% Convolution between the Fourier Transformed
% image and the mask
G = H.*FT_img;
  
% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)  
output_image = real(ifft2(double(G)));
% Displaying Input Image and Output Image
subplot(2, 1, 1), imshow(input_image),
subplot(2, 1, 2), imshow(output_image, [ ]);

% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% IHPF
figure
input_image = imread("img.jpg");  
[M, N] = size(input_image , [1 2]);
FT_img = fft2(im2double(input_image));
D0 = 15; 

% Designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;

% MATLAB library function meshgrid(v, u) returns
% 2D grid which contains the coordinates of vectors
% v and u. Matrix V with each row is a copy 
% of v, and matrix U with each column is a copy of u
[V, U] = meshgrid(v, u);
  
% Calculating Euclidean Distance
D = sqrt(U.^2+V.^2);
  
% Comparing with the cut-off frequency and 
% determining the filtering mask
H = double(D <= D0)-1;
  
% Convolution between the Fourier Transformed
% image and the mask
G = H.*FT_img;
  
% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)  
output_image = real(ifft2(double(G)));
% Displaying Input Image and Output Image
subplot(2, 1, 1), imshow(input_image),
subplot(2, 1, 2), imshow(output_image, [ ]);

% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% BHPF
figure;
input_image = imread("img.jpg");  
[M, N] = size(input_image , [1 2]);
FT_img = fft2(im2double(input_image));
D0 = 15; 
n=2*2;
% Designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;

% MATLAB library function meshgrid(v, u) returns
% 2D grid which contains the coordinates of vectors
% v and u. Matrix V with each row is a copy 
% of v, and matrix U with each column is a copy of u
[V, U] = meshgrid(v, u);
  
% Calculating Euclidean Distance
D = sqrt(U.^2+V.^2);

D = D./ D0;
% Comparing with the cut-off frequency and 
% determining the filtering mask
H = 1./((1+D).^n)-1;
  
% Convolution between the Fourier Transformed
% image and the mask
G = H.*FT_img;
  
% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)  
output_image = real(ifft2(double(G)));
% Displaying Input Image and Output Image
subplot(2, 1, 1), imshow(input_image),
subplot(2, 1, 2), imshow(output_image, [ ]);

% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
originalImage = imread("res.png");   
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

% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% GHPF
figure;
input_image = imread("gs.png");  
[M, N] = size(input_image , [1 2]);
FT_img = fft2(im2double(input_image));
D0 = 15; 
D0 = (D0^2)*2;
% Designing filter
u = 0:(M-1);
idx = find(u>M/2);
u(idx) = u(idx)-M;
v = 0:(N-1);
idy = find(v>N/2);
v(idy) = v(idy)-N;

% MATLAB library function meshgrid(v, u) returns
% 2D grid which contains the coordinates of vectors
% v and u. Matrix V with each row is a copy 
% of v, and matrix U with each column is a copy of u
[V, U] = meshgrid(v, u);
  
% Calculating Euclidean Distance
D = sqrt(U.^2+V.^2);

D = -D.^2;
% Comparing with the cut-off frequency and 
% determining the filtering mask
H = exp(D/D0)-1;
  
% Convolution between the Fourier Transformed
% image and the mask
G = H.*FT_img;
  
% Getting the resultant image by Inverse Fourier Transform
% of the convoluted image using MATLAB library function 
% ifft2 (2D inverse fast fourier transform)  
output_image = real(ifft2(double(G)));
% Displaying Input Image and Output Image
subplot(2, 1, 1), imshow(input_image),
subplot(2, 1, 2), imshow(output_image, [ ]);

% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Laplacian filter
figure
originalImage = imread("contrast.png");  
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

% --- Executes on button press in pushbutton36.
function pushbutton36_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
gray_image = imread("img.jpg");    
subplot(1, 2, 1),
% displaying the image
imshow(gray_image);
title("Original image");

[rows,cols]=size(gray_image);
out=gray_image;
for i=2:rows-1
 for j=2:cols-1
     temp = [gray_image(i-1, j-1) gray_image(i-1, j) gray_image(i-1, j + 1) gray_image(i, j-1) gray_image(i, j) gray_image(i, j + 1) gray_image(i + 1, j-1) gray_image(i + 1, j) gray_image(i + 1, j + 1)];
     temp = sort(temp);
     out(i, j)= temp(5);
end
end
subplot(1, 2, 2),
imshow(out);
title("Median Image");

% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
A=imread('avg.png');
subplot(1,2,1);
imshow(A);
B=rgb2gray(A);

C=double(B);


for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        %Sobel mask for x-direction:
        Gx=((2*C(i+2,j+1)+C(i+2,j)+C(i+2,j+2))-(2*C(i,j+1)+C(i,j)+C(i,j+2)));
        %Sobel mask for y-direction:
        Gy=((2*C(i+1,j+2)+C(i,j+2)+C(i+2,j+2))-(2*C(i+1,j)+C(i,j)+C(i+2,j)));
     
        %The gradient of the image
        %B(i,j)=abs(Gx)+abs(Gy);
        B(i,j)=sqrt(Gx.^2+Gy.^2);
     
    end
end
subplot(1,2,2);
imshow(B); title('Sobel gradient');
