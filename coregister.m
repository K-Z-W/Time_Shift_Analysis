function varargout = V_coregister(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function is used to  register clinical image (e.g. DWI, PWI, ASL_CBF)
% to functional data (mean functional image from ¡®Realign¡¯). Then the TSA
% results could overlay on the clinical image.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Last Modified by WEI WEI, 28-Apr-2016 12:28:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @coregister_OpeningFcn, ...
                   'gui_OutputFcn',  @coregister_OutputFcn, ...
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


function cor_btn_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coregister Button Code  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

 t1=clock;
 handles=guidata(gcf);
 
 sou_filPath=get(handles.sou_edit,'string');
 ref_filPath=get(handles.ref_edit,'string');
 other_filPath=get(handles.other_edit,'string');
 work_dirPath=get(handles.workdir_edit,'string');
 
 filefil_ref=ref_filPath;
 filefil_source=sou_filPath;
 filefil_other=other_filPath;
 

 % \\This part is performing based on SPM\\
 % load coreg_ew;
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {'F:\MRI_packages\TSA\public\home\ccbd019\Desktop\lyt\TIA_data\linux_FunImg\sub011\mean_func.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {'F:\MRI_packages\TSA\public\home\ccbd019\Desktop\lyt\TIA_data\dwi0000.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% batch_coreg = matlabbatch;

 cd(work_dirPath);

 ref = spm_select('ExtList',pwd, ['^' filefil_ref ],inf);
 matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[pwd filesep deblank(ref)]};

 sour = spm_select('ExtList',pwd, ['^' filefil_source ],inf);
 matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[pwd filesep deblank(sour)]};

 if ~isempty(filefil_other) 
     
     other = spm_select('ExtList',pwd, ['^' filefil_other ],inf);
     matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[pwd filesep deblank(other)]};
     
 end
 
 spm_jobman('run',matlabbatch);

% ----------------------------------------------------------------------------
 
 display('Coregister is done');



function coregister_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
guidata(hObject, handles);


function varargout = coregister_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


function workdir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to workdir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function workdir_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function workdir_btn_Callback(hObject, eventdata, handles)

handles=guidata(gcf);
    workdir_Path=uigetdir(pwd);
    if 0~=workdir_Path
        cd(workdir_Path);
        set(handles.workdir_edit,...
            'string',workdir_Path);
    end


function sou_edit_Callback(hObject, eventdata, handles)
% hObject    handle to sou_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function sou_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function sou_dir_btn_Callback(hObject, eventdata, handles)

handles=guidata(gcf);
    sou_fil=uigetfile('*.nii');
    if 0~=sou_fil
        set(handles.sou_edit,...
            'string',sou_fil);
    end


function ref_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ref_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function ref_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ref_dir_btn_Callback(hObject, eventdata, handles)

handles=guidata(gcf)
    ref_fil=uigetfile('*nii');
    if 0~=ref_fil;
        set(handles.ref_edit,...
            'string',ref_fil);
    end


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function edit4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function other_btn_Callback(hObject, eventdata, handles)

handles=guidata(gcf);
    other_fil=uigetfile('*nii');
    if 0~=other_fil
        set(handles.other_edit,...
            'string',other_fil);
    end


function other_edit_Callback(hObject, eventdata, handles)
% hObject    handle to other_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function other_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


