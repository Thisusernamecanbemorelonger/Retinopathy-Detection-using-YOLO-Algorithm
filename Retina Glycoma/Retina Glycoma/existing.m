    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% EDGE DETECTION METHODS %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; %CLC    Clear command window.
clear all; %CLEAR  Clear variables and functions from memory.
close all; %CLOSE Close figure.

%IMREAD Read image from graphics file.
A=imread('C:\Users\vc sir office\Desktop\fv\retina\4.jpg');
imshow(A) %IMSHOW Display image in Handle Graphics figure.  

B = imresize(A, [384 NaN]);%IMRESIZE Resize image.
figure,imshow(B)

C = rgb2gray(B); %RGB2GRAY Convert RGB image or colormap to grayscale.
whos %WHOS List current variables, long form. 
figure,imshow(C)

figure,imhist(C) %IMHIST Display histogram of image data.

%EDGE Find edges in intensity image.
P = edge(C,'prewitt'); %The Prewitt method finds edges using the Prewitt approximation

Q= edge(C,'sobel'); %The Sobel method finds edges using the Sobel approximation 

R= edge(C,'roberts'); %The Roberts method finds edges using the Roberts approximation

S= edge(C,'canny'); %The Canny method finds edges by looking for local maxima of the  gradient 

figure; %FIGURE Create figure window.
subplot (2,2,1); %SUBPLOT Create axes in tiled positions.
imshow(P);
title('Prewitt'); %TITLE  Graph title.

subplot(2,2,2); 
imshow(Q);
title('Sobel');

subplot(2,2,3);
imshow(R);
title('Roberts');

subplot(2,2,4);
imshow(S); 
title('canny'); 

I = double(C);
%CLASS  Return class name of object.
% S = CLASS(OBJ) returns the name of the class of object OBJ.
classType = class(C);
%DOUBLE Convert to double precision.
scalingFactor = double(intmax(classType));
I = I/scalingFactor;

Gx = [-1 1];
Gy = Gx';
%CONV2 Two dimensional convolution.
Ix = conv2(I,Gx,'same');
Iy = conv2(I,Gy,'same');

%NEWFIS Create new FIS.
%   FIS=NEWFIS(FISNAME) creates a new Mamdani-style FIS structure
edgeFIS = newfis('edgeDetection');

%   Add a variable to an FIS.
edgeFIS = addvar(edgeFIS,'input','Ix',[-1 1]);
edgeFIS = addvar(edgeFIS,'input','Iy',[-1 1]);

sx = 0.1; 
sy = 0.1;

%  A membership function can only be added to a variable name
edgeFIS = addmf(edgeFIS,'input',1,'zero','gaussmf',[sx 0]);
edgeFIS = addmf(edgeFIS,'input',2,'zero','gaussmf',[sy 0]);

edgeFIS = addvar(edgeFIS,'output','Iout',[0 1]);

wa = 0.1; wb = 1; wc = 1;
ba = 0; bb = 0; bc = .7;

edgeFIS = addmf(edgeFIS,'output',1,'white','trimf',[wa wb wc]);
edgeFIS = addmf(edgeFIS,'output',1,'black','trimf',[ba bb bc]);

figure
subplot(2,2,1); plotmf(edgeFIS,'input',1); title('Ix');
subplot(2,2,2); plotmf(edgeFIS,'input',2); title('Iy');
subplot(2,2,[3 4]); plotmf(edgeFIS,'output',1); title('Iout')

r1 = 'If Ix is zero and Iy is zero then Iout is white';
r2 = 'If Ix is not zero or Iy is not zero then Iout is black';
r = char(r1,r2);
%PARSRULE Parse fuzzy rules.
edgeFIS = parsrule(edgeFIS,r);
%SHOWRULE Display FIS rules.
showrule(edgeFIS)

Ieval = zeros(size(I));% Preallocate the output matrix
for ii = 1:size(I,1)
    Ieval(ii,:) = evalfis([(Ix(ii,:));(Iy(ii,:));]',edgeFIS);
end

figure; image(Ieval,'CDataMapping','scaled'); colormap('gray');
title('Edge Detection Using Fuzzy Logic')



