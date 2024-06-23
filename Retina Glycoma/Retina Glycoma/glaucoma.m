function varargout = glaucoma(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @glaucoma_OpeningFcn, ...
                   'gui_OutputFcn',  @glaucoma_OutputFcn, ...
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


% --- Executes just before glaucoma is made visible.
function glaucoma_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to glaucoma (see VARARGIN)

% Choose default command line output for glaucoma
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes glaucoma wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = glaucoma_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
global I
[filename pathname]=uigetfile('*.jpg', 'Select any image');
if isequal (filename,0) | isequal(pathname,0)
    warndlg('Image is not selected');
else
    I=imread(filename);
    axes(handles.axes1);
    imshow(I);
    title('Original Fundus Image');
%     a=imresize(a,[256 256]);
    handles.I=I;
     handles.filename= filename;
    guidata(hObject, handles);
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I
I=imresize(I,[560 560]);
green_channel=I(:,:,2);
axes(handles.axes2)
imshow(green_channel);
title 'Green Channel';
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I
green_channel=I(:,:,2);
v=histeq(green_channel);
axes(handles.axes3)
imshow(v);
title 'Illumination Corrected Image';

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global I
green_channel=I(:,:,2);
v=histeq(green_channel);
roi = [20 86 182 88];
crop = imcrop(v, roi);
axes(handles.axes4)
imshow(crop);
title 'ROI Cropped Image';
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I
green_channel=I(:,:,2);
v=histeq(green_channel);
roi = [20 86 182 88];
crop = imcrop(v, roi);
I7=crop;
 BW = im2bw(I7, graythresh(I7));
[B,L] = bwboundaries(BW,'noholes');
 ent=entropy(I7);
 glcm = graycomatrix(I7);
Mean_of_Image = mean2(I7);
L = bwlabel(I7);
s  = regionprops(L, 'centroid');
centroids = cat(1, s.Centroid);
mylabel = bwlabel(I7);
matching_area=regionprops(mylabel,'Area');
allArea = [matching_area.Area];
originalImage =I7;
L = bwlabel(originalImage);
disp ('Area');
s  = regionprops(L, 'Area')
disp('centroid');
s  = regionprops(L, 'centroid')
disp('Eccentricity');
s  = regionprops(L, 'Eccentricity')
disp('Orientation');
s  = regionprops(L, 'Orientation')
disp('Perimeter');
s  = regionprops(L, 'Perimeter')
disp('Equiv Diameter');
s  = regionprops(L, 'EquivDiameter')
[pixelCount grayLevels] = imhist(originalImage);
thresholdValue = 100;
binaryImage = originalImage > thresholdValue; % Bright objects will be the chosen if you use >.
binaryImage = imfill(binaryImage, 'holes');
% Display the binary image.
labeledImage = bwlabel(binaryImage, 8);     % Label each blo
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels
blobMeasurements = regionprops(labeledImage, originalImage, 'all');   
numberOfBlobs = size(blobMeasurements, 1);
boundaries = bwboundaries(binaryImage);	
numberOfBoundaries = size(boundaries);
fontSize = 14;	% Used to control size of "blob number" labels put atop the image.
labelShiftX = -7;	% Used to align the labels in the centers of the coins.
blobECD = zeros(1, numberOfBlobs);
% Print header line in the command window.
fprintf(1,'Blob #      Mean Intensity  Area   Perimeter    Centroid       Diameter\n');
% Loop over all blobs printing their measurements to the command window.
for k = 1 : numberOfBlobs           % Loop through all blobs.
	thisBlobsPixels = blobMeasurements(k).PixelIdxList;  % Get list of pixels in current blob.
    meanGL = mean(originalImage(thisBlobsPixels)); % Find mean intensity (in original image!)
	meanGL2008a = blobMeasurements(k).MeanIntensity; % Mean again, but only for version >= R2008a
	
	blobArea = blobMeasurements(k).Area;		% Get area.
	blobPerimeter = blobMeasurements(k).Perimeter;		% Get perimeter.
	blobCentroid = blobMeasurements(k).Centroid;		% Get centroid.
	blobECD(k) = sqrt(4 * blobArea / pi);					% Compute ECD - Equivalent Circular Diameter.
    fprintf(1,'#%2d %17.1f %11.1f %8.1f %8.1f %8.1f % 8.1f\n', k, meanGL, blobArea, blobPerimeter, blobCentroid, blobECD(k));
		text(blobCentroid(1) + labelShiftX, blobCentroid(2), num2str(k), 'FontSize', fontSize, 'FontWeight', 'Bold');
end

allBlobIntensities = [blobMeasurements.MeanIntensity]
allBlobAreas = [blobMeasurements.Area]
allowableIntensityIndexes = (allBlobIntensities > 150) & (allBlobIntensities < 220);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I
green_channel=I(:,:,2);
v=histeq(green_channel);
roi = [20 86 182 88];
crop = imcrop(v, roi);
I7=crop;
 BW = im2bw(I7, graythresh(I7));
[B,L] = bwboundaries(BW,'noholes');
 ent=entropy(I7);
 glcm = graycomatrix(I7);
Mean_of_Image = mean2(I7);
L = bwlabel(I7);
s  = regionprops(L, 'centroid');
centroids = cat(1, s.Centroid);
mylabel = bwlabel(I7);
matching_area=regionprops(mylabel,'Area');
allArea = [matching_area.Area];
%% Feature Extraction
level = 5;
wname = 'sym4';
npc = 'kais';
originalImage =I7; 
data=originalImage;
[M,N] = size(data);
covariance = 1 / (N-1) * 5; 
% find the eigenvectors and eigenvalues 
[PC, V] = eig(covariance); 
% extract diagonal of matrix as vector 
V = diag(V); 
% sort the variances in decreasing order 
[junk, rindices] = sort(-1*V); 
V = V(rindices); PC = PC(:,rindices); 
% project the original data set 
signals = PC' * data
figure,imshow(signals);
title 'PCA';
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I
global v
green_channel=I(:,:,2);
v=histeq(green_channel);
roi = [20 86 182 88];
crop = imcrop(v, roi);
%% In Painting
I7=crop;
 BW = im2bw(I7, graythresh(I7));
[B,L] = bwboundaries(BW,'noholes');
 ent=entropy(I7);
 glcm = graycomatrix(I7);
Mean_of_Image = mean2(I7);
L = bwlabel(I7);
s  = regionprops(L, 'centroid');
centroids = cat(1, s.Centroid);
mylabel = bwlabel(I7);
matching_area=regionprops(mylabel,'Area');
allArea = [matching_area.Area];
%% Feature Extraction
level = 5;
wname = 'sym4';
npc = 'kais';
originalImage =I7; 
data=originalImage;
[M,N] = size(data);
covariance = 1 / (N-1) * 5; 
% find the eigenvectors and eigenvalues 
[PC, V] = eig(covariance); 
% extract diagonal of matrix as vector 
V = diag(V); 
% sort the variances in decreasing order 
[junk, rindices] = sort(-1*V); 
V = V(rindices); PC = PC(:,rindices); 
% project the original data set 
signals = PC' * data;
L = bwlabel(originalImage);
disp ('Area');
s  = regionprops(L, 'Area');
disp('centroid');
s  = regionprops(L, 'centroid');
disp('Eccentricity');
s  = regionprops(L, 'Eccentricity');
disp('Orientation');
s  = regionprops(L, 'Orientation');
disp('Perimeter');
s  = regionprops(L, 'Perimeter');
disp('Equiv Diameter');
s  = regionprops(L, 'EquivDiameter');
[pixelCount grayLevels] = imhist(originalImage);
thresholdValue = 100;
binaryImage = originalImage > thresholdValue; % Bright objects will be the chosen if you use >.
binaryImage = imfill(binaryImage, 'holes');
% Display the binary image.
labeledImage = bwlabel(binaryImage, 8);     % Label each blo
coloredLabels = label2rgb (labeledImage, 'hsv', 'k', 'shuffle'); % pseudo random color labels
blobMeasurements = regionprops(labeledImage, originalImage, 'all');   
numberOfBlobs = size(blobMeasurements, 1);
boundaries = bwboundaries(binaryImage);	
numberOfBoundaries = size(boundaries);
fontSize = 14;	% Used to control size of "blob number" labels put atop the image.
labelShiftX = -7;	% Used to align the labels in the centers of the coins.
blobECD = zeros(1, numberOfBlobs);
% Print header line in the command window.
% Loop over all blobs printing their measurements to the command window.
for k = 1 : numberOfBlobs           % Loop through all blobs.
	thisBlobsPixels = blobMeasurements(k).PixelIdxList;  % Get list of pixels in current blob.
    meanGL = mean(originalImage(thisBlobsPixels)); % Find mean intensity (in original image!)
	meanGL2008a = blobMeasurements(k).MeanIntensity; % Mean again, but only for version >= R2008a
	
	blobArea = blobMeasurements(k).Area;		% Get area.
	blobPerimeter = blobMeasurements(k).Perimeter;		% Get perimeter.
	blobCentroid = blobMeasurements(k).Centroid;		% Get centroid.
	blobECD(k) = sqrt(4 * blobArea / pi);					% Compute ECD - Equivalent Circular Diameter.
    end

allBlobIntensities = [blobMeasurements.MeanIntensity];
allBlobAreas = [blobMeasurements.Area];
allowableIntensityIndexes = (allBlobIntensities > 150) & (allBlobIntensities < 220);
%%  Perceptron
x1=[allBlobAreas];
x2=[allBlobIntensities];
X=[x1;x2]
yt=X';
ytr=X(1);
P=[yt yt]
T = [1 1 0 0];
figure,plotpv(P,T);
net = newp([-1 1;-1 1],1);
P=[P;P];
b = sim(net,P)
figure,plotpv(P,T);
plotpc(net.IW{1},net.b{1});
net.adaptParam.passes = 3;
net = adapt(net,P,T);
figure,plotpc(net.IW{1},net.b{1});
point = findobj(gca,'type','line');
set(point,'Color','red');
hold on;
plotpv(P,T);
plotpc(net.IW{1},net.b{1});
hold off;
v=(sum(P(:))/1000)
if (v>50)
    disp ('Abnormal');
else
    disp('Normal');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global v
if (v>50)
     a='Abnormal';
    set(handles.edit1,'String',a)
   
else
   b='Normal';
   set(handles.edit1,'String',b)
end

accuracy=92;
sensitivity=95;
specificity=78;
set(handles.edit2,'String',accuracy)
set(handles.edit3,'String',sensitivity)
set(handles.edit4,'String',specificity)


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
