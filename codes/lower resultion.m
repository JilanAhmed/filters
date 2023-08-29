A = imread("img.jpg");
subplot (1,2,1);
imshow(A)
sz = size(A);
xg = 1:sz(1);
yg = 1:sz(2);
F = griddedInterpolant({xg,yg},double(A));xq = (0:5/6:sz(1))';
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