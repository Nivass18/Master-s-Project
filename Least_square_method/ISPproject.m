function varargout = ISPproject(varargin)
% ISPPROJECT MATLAB code for ISPproject.fig
%      ISPPROJECT, by itself, creates a new ISPPROJECT or raises the existing
%      singleton*.
%
%      H = ISPPROJECT returns the handle to a new ISPPROJECT or the handle to
%      the existing singleton*.
%
%      ISPPROJECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK ISPproject ISPPROJECT.M with the given input arguments.
%
%      ISPPROJECT('Property','Value',...) creates a new ISPPROJECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ISPproject_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ISPproject_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ISPproject

% Last Modified by GUIDE v2.5 25-Jan-2019 00:34:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ISPproject_OpeningFcn, ...
                   'gui_OutputFcn',  @ISPproject_OutputFcn, ...
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
end

% --- Executes just before ISPproject is made visible.
function ISPproject_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined ISPproject a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ISPproject (see VARARGIN)

% Choose default command line output for ISPproject
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ISPproject wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = ISPproject_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined ISPproject a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined ISPproject a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

n = str2num(get(hObject,'String'));
setappdata(0,'n',n);

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined ISPproject a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press ISPproject main.
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined ISPproject a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

N = str2num(get(hObject,'String'));
setappdata(0,'N',N);

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
end


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined ISPproject a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function equation_Callback(hObject, eventdata, handles)
% hObject    handle to equation (see GCBO)
% eventdata  reserved - to be defined ISPproject a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA
% Hints: get(hObject,'String') returns contents of equation as text
%str2double(get(hObject,'String')) returns contents of equation as a double
end

% --- Executes during object creation, after setting all properties.
function equation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to equation (see GCBO)
% eventdata  reserved - to be defined ISPproject a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function main_Callback(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  reserved - to be defined ISPproject a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

n = getappdata(0,'n');% The number set of data that
...we want to analysis the lsq for the graph
    %n= str2num(get(handle.edit2,'string'));

N = getappdata(0,'N');%The number of poly fit that you want
%N= str2num(get(hadle.edit2,'string'));
if (isempty(n)==0)&&(isempty(N)==0)
   if (isstring(n)==0)&&(isstring(N))&&(n>0)&&(N>0)
       if(isfloat(n)==0)&&(isfloat(N)==0)
        

filename=uigetfile({'*.csv'},'File selector');%Reading the x and y values
Q=dlmread(filename);x=Q(1:n,1);y=Q(1:n,2);
fprintf('Reading x value\nReading y value\n');

% Calculation of S-Vector
%FOR THE EXPLANATION OF THE BELOW FORLOOP REFER-1.1 IN DOCUMENT
for i = 1:(2*N+1)
    s=0;
    for j=1:n
        s = s+(x(j))^(i-1);S(i)=s;
    end
end

%Calculation of T vector
%FOR THE EXPLANATION OF THE BELOW FORLOOP REFER-1.2 IN DOCUMENT
for i=1:(N+1)
    t=0;
    for j=1:n
        t= t + (y(j)*(x(j))^(i-1)); T(i)= t;
    end
end      

%%Matrics arrangement%%
%FOR THE EXPLANATION OF THE BELOW FORLOOP REFER-1.3 IN DOCUMENT
A=zeros(N+1,N+1);
k=0;
for i=1:(N+1)
        k=i-1;
    for j=1:(N+1)
        a=k+2;
        switch a
            case (i+j)
                A(i,j)=S(k+1);          
        end
        if(k<2*N)
           k=k+1;
        end
    end
end

%FOR THE EXPLANATION OF THE BELOW FORLOOP REFER-1.4 IN DOCUMENT
%Arrangement of T-Vecot in Matrics B
for i=1:(N+1)
    B(i)=T(i);
end

%Calculation of co-efficient of the polynomial equation.
X=(A^-1)*B';

%The equation is printed
fprintf('the Equation is :\n');
D=' ';
for i=1:(N+1)
    if((i-1)==0)
        fprintf('%+d',X(i));
    else
        fprintf('%+d*x^%d',X(i),(i-1));
    end
    if(X(i)>0)
        L='+';
    else
        L='';
    end
    U=num2str(X(i))*str2sym('x').^num2str((i-1));
    M=char(U);
    %concatenation of the equation for printing it in the gui window
     D=strcat(L,num2str(M),D);
end
fprintf('\n');
set(handles.equation,'String',D); %setting equation as the output data

%displaying the x and y values taken as input
fprintf('x-value\t\ty-values\n');
for i=1:n
    fprintf('\t%d',x(i));
    fprintf('\t\t\t\t%d',y(i));
    fprintf('\n')
end

%The A matrics is printed
fprintf('the matrixs of A:\n');
for i=1:(N+1)
    for j=1:(N+1)
        fprintf('%d,',A(i,j));
    end
     fprintf('\n');
end

%printing the B and X values
fprintf('The Matrics B\n');fprintf('%d\t',B);fprintf('\n');
fprintf('The Matrics X\n');fprintf('%d\t',X);fprintf('\n');

%for loop for calculating the polyfit-curves points.
for i=1:n
    V=0;
    for j=1:(N+1)
        V=V+(X(j)*(x(i)^(j-1)));
    end
    R(i)=V; W(i)=x(i);
end

%syntaxs for plotting the curve
plot(W,R,x,y,'o');
xlabel('X-values');
ylabel('Y-Values');
        
    else
        errordlg('you have not entered a proper value ');
       end
 else
    errordlg('you have not entered a proper value');
   end
  else
    errordlg('please enter the value for empty block or you must have entered a string');
end
end
