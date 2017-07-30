function varargout =automask(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to create a brain mask for your own data which will
% be used in the time shift analysis.
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Citing Information:
% Automask performs based on w_Automask and w_Cliplevel, that were written
% by Xin-di Wang.

% State Key Laboratory of Cognitive Neuroscience and Learnling,
% Beijing Normal University, Beijing, PR China
% sandywang.rest@gmail.com
%                                    
% Modified according to AFNI's 3dAutomask.

%--------------------------------------------------------------------------
% Last Modified by WEI WEI, 29-Apr-2016 10:12:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @automask_OpeningFcn, ...
                   'gui_OutputFcn',  @automask_OutputFcn, ...
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


function automask_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

guidata(hObject, handles);


function varargout = automask_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


function mean_dir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mean_dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function mean_dir_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function mask_dir_btn_Callback(hObject, eventdata, handles)

handles=guidata(gcf);
    
%     DPATH=uigetdir(pwd); 
% 
%     if 0~=DPATH
%         cd(DPATH);
%         DirectoryPath=uigetfile( {'*.nii';'*.img';'*.hdr'},'File Selector');
%         workPath=[DPATH filesep DirectoryPath ];
%         set(handles.mean_dir_edit,...
%             'string',workPath);
%     end

[filename,pathname]= uigetfile(...
                                                       {'*.nii;*.img'});
    str_A=strcat(pathname,filename);
    if 0 ~= str_A
        cd(pathname);
        set(handles.mean_dir_edit,...
              'string',str_A);
    end


function AM_B_Callback(hObject, eventdata, handles)

 handles=guidata(gcf);
 workPath=get(handles.mean_dir_edit,'string');
 
 w_Automask(workPath,'fun_mask',1);
 
 display('AutoMask is done !');
