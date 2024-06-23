function varargout =cellcountGUI2(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cellcountGUI2_OpeningFcn, ...
                   'gui_OutputFcn',  @cellcountGUI2_OutputFcn, ...
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
% End initialization code


% --- Executes just before cellcountGUI2 is made visible.
function cellcountGUI2_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cellcountGUI2 (see VARARGIN)

% Choose default command line output for cellcountGUI2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cellcountGUI2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cellcountGUI2_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse.

function browse_Callback(~, ~, handles)
global im;
global imOld;
global map;
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname]=uigetfile({'*.tif';'*.png';'*.tiff';'*.jpg';'*.bmp';'*.gif'},'Open File');
[im,map]=imread(filename);
imOld=im;
axes(handles.axes1);
imshow(im,map);
axis off
title(filename,'fontsize',16)
set(handles.edit1,'string','');

% --- Executes on button press in count.

function count_Callback(~, ~, handles)
% hObject    handle to count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im;
global imOld;
global imOut;
global level;
global overlay1;
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~ , ~, c]=size(im);
if(c>1) % color image
switch get(handles.popupmenu1,'Value')
    case 2 % Red
        im1=im(:,:,1);
    case 3 % Green
        im1=im(:,:,2);
    case 4 % Blue
        im1=im(:,:,3);
     otherwise
        im1=im(:,:,2);
end
else % grayscale image
    im1=im;
end
% frequency filtering (band pass)
[row, column]=size(im1);
imFFT=fft2(im1,row,column);  % Fourier transformed image, same size as original image
imFFT=fftshift(imFFT);  % Fourier transformed image centered
temp=min(row,column);

temp1=max(row,column);
% Build a bandpass filter.  Everything within the pass band is one,
% everything outside is 0.0.
ri=repmat([1:temp]',1,temp); % the size of the original image is 1200X1600
ci=repmat([1:temp],temp,1);
ri=ri-temp/2.0;
ci=ci-temp/2.0;
r(1:temp,(temp1-temp)/2+1:(temp1+temp)/2)=sqrt(ri.*ri+ci.*ci);
r(1:temp,1:(temp1-temp)/2)=max(r(:)).*ones(temp,(temp1-temp)/2); % pad the margin
r(1:temp,(temp1+temp)/2+1:temp1)=r(1:temp,1:(temp1-temp)/2);
radius1=6;
radius2=30;
index= r<radius2 & r>radius1;
filt=zeros(row,column); % everything outside the pass band is 0
filt(index)=1.0;
% Apply the filter
imFFT2=imFFT.*filt; % Masked centered fft of image with bandpass filter
imOut=ifft2(fftshift(imFFT2));
imOut=real(imOut);
axes(handles.axes2);
if(min(imOut(:))<0) % make all values positive
    imOut=imOut+abs(min(imOut(:)));
end
imOut=imOut./max(imOut(:));
level=graythresh(imOut);
BW=im2bw(imOut,level);
bw2 = imfill(BW,'holes');
bw3 = imopen(bw2, ones(5,5));
bw4 = bwareaopen(bw3, 40);
bw4_perim = bwperim(bw4,8);
SE=[0 1 0
    1 1 1
    0 1 0];
bw4_perim=imdilate(bw4_perim,SE);
colorV=get(handles.popupmenu1,'Value');
if(colorV==2)
    rgbV=[0.3 1 0.3];
else
    rgbV=[1 0.3 0.3];
end
overlay1=imoverlay(imOld,bw4_perim,rgbV);
axes(handles.axes2)
imshow(overlay1)                    
[~,n1]=bwlabel(BW);
title(['N=',num2str(n1)],'color','red','fontsize',16)


% --- Executes on slider movement.
function slider1_Callback(~, ~, handles)
global imOld;
global imOut;
global level;
global overlay1;
global BW;
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of
%        slider
threshold=get(handles.slider1,'value');
set(handles.edit1,'string',num2str(threshold));
BW=im2bw(imOut,level+threshold);
bw2 = imfill(BW,'holes');
bw3 = imopen(bw2, ones(5,5));
bw4 = bwareaopen(bw3, 40);
bw4_perim = bwperim(bw4,8);
SE=[0 1 0
    1 1 1
    0 1 0];
bw4_perim=imdilate(bw4_perim,SE);
colorV=get(handles.popupmenu1,'Value');
if(colorV==2)
    rgbV=[0.3 1 0.3];
else
    rgbV=[1 0.3 0.3];
end
overlay1=imoverlay(imOld,bw4_perim,rgbV);
axes(handles.axes2)
imshow(overlay1)                    
[~,n1]=bwlabel(BW);
title(['N=',num2str(n1)],'color','red','fontsize',16)
drawnow


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, ~, ~)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function edit1_Callback(~, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, ~, ~)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in output.
function output_Callback(~, ~, ~)
global overlay1;
global map;
global BW;
% hObject    handle to output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure
imshow(overlay1,map)                    
[~,n1]=bwlabel(BW);
title(['N=',num2str(n1)],'color','red','fontsize',16)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(~, ~, ~)
global im;
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
se=strel('disk',20);
im=imtophat(im,se);

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, ~, ~)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

