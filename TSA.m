function varargout = TSA(varargin)

% Toolbox for Time Shift Analysis of fMRI data
% Release = Version 1.1
% Copyright(c)
% Written by Wei Wei on May-2016
% Zhejiang Key Laboratory for Research in Assessment of Cognitive Impairments,Center for Cognition and Brain Disorder, Hangzhou Normal University, China
% Mail to author: <"weilisten@163.com"> Wei Wei

%--------------------------------------------------------------------------

% Citing Information:
% The processing was carried out by using REST-Time Shift Analysis, which is based on SPM12 and RESTplus.

% Last Modified by WEI WEI, 24-July-2017 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TSA_OpeningFcn, ...
                   'gui_OutputFcn',  @TSA_OutputFcn, ...
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



function TSA_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = TSA_OutputFcn(hObject, eventdata, handles) 


% Get default command line output from handles structure
varargout{1} = handles.output;



function Directory_btn_Callback(hObject, eventdata, handles)

handles=guidata(gcf);
    DirectoryPath=uigetdir(pwd);
    if 0~=DirectoryPath
        cd(DirectoryPath);
        set(handles.Directory_edit,...
            'string',DirectoryPath);
    end



function slice_ord_Callback(hObject, eventdata, handles)

 SliceOrder=get(hObject,'string');
     SliceOrderText =eval(['[',SliceOrder,']']);
     set(hObject,'string',num2str(SliceOrderText))
 


% --- Executes on button press in pre_btn.
function pre_btn_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Preprocessing  Module     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 t1=clock;
 handles=guidata(gcf);
 DirectoryPath=get(handles.Directory_edit,'string');
 Workdirpath=inpath_Misc(DirectoryPath,'GetParentPath');
 
 mkdir( Workdirpath,'CovariatesParameter');
 OutDirCovPath=[Workdirpath filesep 'CovariatesParameter'];
 

 CutNumber=10;
 TRstr=get(handles.TR_edit,'string');
 TR=eval(['[',TRstr,']']);
 Slice_num_str=get(handles.slice_num_edit,'string');
 SliceNumber=eval(['[',Slice_num_str,']']);
 Slice_Ord_str=get(handles.slice_ord,'string');
 SliceOrder=eval(['[',Slice_Ord_str,']']);
 Ref_Slice_str=get(handles.ref_edit,'string');
 ReferenceSlice=eval(['[',Ref_Slice_str,']']);
 
 IsDicomToNifti=get(handles.d2n_checkbox,'value');
 IsSliceTiming=get(handles.ST_checkbox,'value');

 
 RTP_num_str=get(handles.RTP_edit,'string');
 RTPnumber=eval(['[',RTP_num_str,']']);
 Low_str=get(handles.low_edit,'string');
 Low=eval(['[',Low_str,']']);
 High_str=get(handles.high_edit,'string');
 High=eval(['[',High_str,']']);
 Smooth_str=get(handles.smooth_edit,'string');
 FWHM=eval(['[',Smooth_str,']']);
 

 %   \\ Preprocessing Checkboxes \\
 
 if 1==IsDicomToNifti 
     if 1==IsSliceTiming
         OperationList={'EPIDICOMTONIFTI','RemoveFirstTimePoints','SliceTiming','Realign','Smooth','Detrend','RegressOutCovariates','Filter'};
         Pf=get_Postfix('DicomToNifti'); 
         Pd=get_Postfix('SliceTiming');
     else
         OperationList={'EPIDICOMTONIFTI','RemoveFirstTimePoints','Realign','Smooth','Detrend','RegressOutCovariates','Filter'};
         Pf=get_Postfix('DicomToNifti'); 
         Pd='';
     end
     
 else
     if 1==IsSliceTiming
     OperationList={'RemoveFirstTimePoints','SliceTiming','Realign','Smooth','Detrend','RegressOutCovariates','Filter'};
     Pf=''; 
     Pd=get_Postfix('SliceTiming');
     else
     OperationList={'RemoveFirstTimePoints','Realign','Smooth','Detrend','RegressOutCovariates','Filter'};
     Pf=''; 
     Pd='';
     end
     
 end
 
% ----------------------------------------------------------------------------------- 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MAIN    PREPROCESSING    PART            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% EpiDicomtoNifti
    InputParameter.EpiDicomToNifti.InDirFunRaw=DirectoryPath;
    InputParameter.EpiDicomToNifti.OutDirFunImg=[DirectoryPath Pf];
    
% pipeline_batch({'EPIDICOMTONIFTI'},InputParameter)
%%
% Remove Time Points

  InputParameter.RemoveFirstTimePoints.InDirFunImg=InputParameter.EpiDicomToNifti.OutDirFunImg;
  InputParameter.RemoveFirstTimePoints.OutDirFunImg=[InputParameter.RemoveFirstTimePoints.InDirFunImg,...
                                                       get_Postfix('RemoveFirstTimePoints')];
  InputParameter.RemoveFirstTimePoints.TimePointsAmount=RTPnumber;

% pipeline_batch({'RemoveFirstTimePoints'},InputParameter)
%%
% Slice Timing
 
 InputParameter.SliceTiming.InDirFunImg=InputParameter.RemoveFirstTimePoints.OutDirFunImg;
 InputParameter.SliceTiming.OutDirFunImg=[InputParameter.SliceTiming.InDirFunImg,...
                                           Pd];
 
 InputParameter.SliceTiming.SliceNumber=SliceNumber;
 InputParameter.SliceTiming.SliceOrder=SliceOrder;
 InputParameter.SliceTiming.ReferenceSlice=ReferenceSlice;
 InputParameter.SliceTiming.TR=TR;
 
% pipeline_batch({'SliceTiming'},InputParameter);
%%
% Realign
 
 InputParameter.Realign.InDirFunImg= InputParameter.SliceTiming.OutDirFunImg;
 InputParameter.Realign.OutDirFunImg=[InputParameter.Realign.InDirFunImg,...
                                       get_Postfix('Realign')];
 InputParameter.Realign.RealignParameterDir=[Workdirpath filesep get_Postfix('RealignParameter')];
 
% pipeline_batch({'Realign'},InputParameter);
 
%%
 % Smooth
 
 InputParameter.Smooth.Gaussian.InDirFunImg=InputParameter.Realign.OutDirFunImg;
 InputParameter.Smooth.Gaussian.OutDirFunImg=[InputParameter.Smooth.Gaussian.InDirFunImg,...
                                             get_Postfix('Smooth')];
 InputParameter.Smooth.Gaussian.FWHM=FWHM;    

% pipeline_batch({'Smooth'},InputParameter);

%%
% Detrend

 InputParameter.Detrend.InDirFunImg=InputParameter.Smooth.Gaussian.OutDirFunImg;
 InputParameter.Detrend.OutDirFunImg=[InputParameter.Detrend.InDirFunImg,...
                                    get_Postfix('Detrend')];
 InputParameter.Detrend.CutNumber=CutNumber;
    
% pipeline_batch({'Detrend'},InputParameter);

%%
% RegressOutCovariates
 InputParameter.RegressOutCovariates.InDirFunImg=InputParameter.Detrend.OutDirFunImg;
 InputParameter.RegressOutCovariates.OutDirFunImg=[InputParameter.RegressOutCovariates.InDirFunImg get_Postfix('RegressOutCov')];
 InputParameter.RegressOutCovariates.OutDirCov= OutDirCovPath;
 
 InputParameter.RegressOutCovariates.IsRemoveIntercept= 0;
 InputParameter.RegressOutCovariates.PolynomialTrend=1;
 InputParameter.RegressOutCovariates.IsWholeBrain= 0;
 InputParameter.RegressOutCovariates.IsCSF = 0;
 InputParameter.RegressOutCovariates.IsWhiteMatter = 0;
 InputParameter.RegressOutCovariates.IsHeadMotion_Rigidbody6= 1;
 InputParameter.RegressOutCovariates.IsOtherCovariatesROI = 0;
 InputParameter.RegressOutCovariates.InDirRealignParameter =InputParameter.Realign.RealignParameterDir ;
 InputParameter.RegressOutCovariates.OtherCovariatesROIList = '';

% pipeline_batch({'RegressOutCovariates'},InputParameter);

%%
% Filter
 InputParameter.Filter.InDirFunImg=InputParameter.RegressOutCovariates.OutDirFunImg;
 InputParameter.Filter.OutDirFunImg=[InputParameter.Filter.InDirFunImg get_Postfix('Filter')];
 InputParameter.Filter.InFileMask='';
 InputParameter.Filter.SamplePeriod=TR;
 InputParameter.Filter.LowPass_HighCutoff=High;
 InputParameter.Filter.HighPass_LowCutoff=Low;
 InputParameter.Filter.IsAddMeanBack='Yes';
 InputParameter.Filter.CutNumber=CutNumber;  

% pipeline_batch({'Filter'},InputParameter);
 
%%
% OperationList
 pipeline_batch(OperationList,InputParameter);

 
%% Auto Mask

 mkdir(Workdirpath,'mask');
 mask_dir=[Workdirpath filesep 'mask'];
 RP_dir=[ Workdirpath filesep 'RealignParameter'];
 mkdir(Workdirpath,'HeadMotion');
 MotionDir=[Workdirpath filesep 'HeadMotion'];
 
 cd(RP_dir)
%  sub_dir=dir('sub*');
Exclude_file=dir('Ex*');
headmotion_file=dir('Head*');

movefile( Exclude_file.name , MotionDir );
movefile( headmotion_file.name , MotionDir );


sub_dir=dir;
 
for h=3:length(sub_dir)
     cd(mask_dir)
     mkdir(sub_dir(h).name);
     mean_dir=[RP_dir filesep sub_dir(h).name];
     cd(mean_dir)
     mean_file=dir('mean*');

     copy_dir=[mask_dir filesep sub_dir(h).name];
     copyfile(mean_file.name, copy_dir);
 
     cd(copy_dir)
     old=dir('mean*');
     new_name='mean_origin.nii';
 
     movefile(old.name, new_name);
   
     % 3dAutomask

     f_mask=[copy_dir,'\mean_origin.nii'];
     w_Automask(f_mask,'fun_mask',1);
 end
 display('Automask is done');
 
 display('Preprocessing is done');
 
 


   
function tsa_btn_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Main TSA Module     %
%%%%%%%%%%%%%%%%%%%%%%%%%
 

%% Define work directory

handles=guidata(gcf);
 DirectoryPath=get(handles.Directory_edit,'string');
 Workdirpath=inpath_Misc(DirectoryPath,'GetParentPath');

 IsTSA=get(handles.tsa_checkbox,'value');
 IsPreprocessing=get(handles.pre_checkbox,'value');
 
 if 1==IsPreprocessing             % Perform TSA after completing preprocessing with this tool
    if 1==IsTSA                    
        cd(Workdirpath)
        dataDir=dir('*RSDCF');
     
        DataDir=[Workdirpath filesep dataDir.name];
        
    end
    
 else                              % Perform TSA with preprocessed fMRI data 
    if 1==IsTSA                    
        
        DataDir=DirectoryPath;
 
    end
end


%% Calculate correlation between shifted time course of the voxel and nean time course
 
 SRstr=get(handles.SR_edit,'string');
 SR=eval(['[',SRstr,']']);
 HMstr=get(handles.HM_edit,'string');
 
  TRstr=get(handles.TR_edit,'string');
  TR=eval(['[',TRstr,']']);
 
 
    datadir =  DataDir;            % the path of the preprocessed resting-state data
%     dirname = 'sub*';              % the prefix of the subjects' directory
%     filename = 'sub';              % the prefix of the preprocessed resting-state data, this prefix should be same for all subjects

    maskfile_dir=[Workdirpath filesep 'mask'];
    cd(maskfile_dir)
%     sub_dir=dir('sub*');
    sub_dir=dir;

for l=3:length(sub_dir)
    maskfilePath=[maskfile_dir filesep sub_dir(l).name];
    maskfile=[maskfilePath filesep 'fun_mask.nii'];


    shift_range = SR;    % set the time shift range

    cd(datadir)
%   dirlist=dir(dirname);   
    dirlist=dir;     % list of subjects
    dirlist_l=length(dirlist);
    
    mkdir(datadir,'time_shift_results');          % make the directory of output file
    outputdir=[datadir filesep 'time_shift_results'];     % the directory name of output file


if dirlist_l==2
    disp('There is no subject in the working directory!')
else

        cd (outputdir)
        mkdir (dirlist(l).name);
        cd (dirlist(l).name)
        destidir=pwd;
        
        outputdir2=[datadir filesep dirlist(l).name];
        cd (outputdir2)
        
       % \\Select the preprocessed resting-state data\\
         File_name = dir;
       % ----------------------------------------------
         
         if length(File_name)==2
            disp('There is no file in the subject directory!')
         else
                 
        % \\Read the time courses of all brain voxels (origin_voxel)\\
           
           V = spm_vol(File_name(3).name);                      % read the header information (# num of fname is num of volume #)
           vn=size(V,1);                                % volume no
           
           origin_voxel= V_get_tc_roi(V, maskfile);     % Read the time course of whole brain mask
        % -----------------------------------------------------------
        
        
        % \\ If you want to use own reference mask \\
            
             if length(HMstr)~=7
                hemi_voxel=V_get_tc_roi(V,HMstr);
             else                 
                hemi_voxel=origin_voxel;
             end
         % ----------------------------------------
             
             hemi_ave=mean(hemi_voxel,2);              % Mean of the brain time courses
             
             [A B C]=rp_readfile(maskfile);            % A:data file; B:size of voxel; C:header information
             [x,y,z]=ind2sub(size(A),find(A~=0));
             hemi_ave_ml=repmat(hemi_ave,1,size(x,1));

             if isempty(x)
                 disp('the mask has no nonzero voxel!')
             else
                 copyfile (maskfile,destidir);
                 m=shift_range;
                     for s=-m:m
                         D=zeros(size(A));
                         R = V_Cor_Pair(hemi_ave_ml(1+m:vn-m,:)',origin_voxel((1+m-s):(vn-m-s),:)'); % calculate the correlation use the function 'V_Cor_Pair'
                         for i=1:size(x,1)
                             D(x(i),y(i),z(i))=R(i);
                         end
                                                  
                          str=['M',num2str(s+m+1)];
                          eval([str,'=D;']);
                          eval([str,'=',str,'(:);']);
                          eval([str,'(find(', str,'==0))=[];']);
                                     
                         clear D;
                     end
                     
                      M = [];
                          for i = 1:s+m+1
                              stri=num2str(i);
                              eval(['M=[M,M',stri,'];']);
                          end
                      
                       M=M';
                       
              
 %% Find the maximal correlation value,put the time shift value which showed maximal corrlation to the voxel
 
 
                       cd([outputdir filesep dirlist(l).name])
                       conn_voxel=M;
                       
                       D=zeros(size(A));        % D is the matrix to save the time shift value showes maximal correlation
                       E=zeros(size(A));        % E is the contrary matrix of D.
                                                          % E is used for display (the time shift values in lesion areas are smaller than zero, make them positive is easy to show the results)
                       
                       
                       m=shift_range;
                       
                           for i=1:size(x,1)
                            conn_voxel(isnan(conn_voxel)==1)=0;
                            if max(conn_voxel(:,i))==0
                                tmp(i)=m+1;
                            else
                                tmp(i)=find(conn_voxel(:,i)==max(conn_voxel(:,i)));
                            end
                            D(x(i),y(i),z(i))=(tmp(i)-m-1)*TR;
                           end
                           
                           E=-D;
                           
     % Save the time shift value which shows maximal and write the result matrix into nifti file
     
                        newname1=strcat(num2str(m),'Wbhemi_ts');
                        rp_writefile(D,newname1,size(A),B,C,'double');
                        newname2=strcat(num2str(m),'Wbhemi_ts_oppo');
                        rp_writefile(E,newname2,size(A),B,C,'double');
                        
                        clear D E
                        
                        display([dirlist(l).name blanks(1) 'is done'])
                        
                       

             end
         end

end


end

 % Move TSA results to output directory
 
    cd(DataDir)

    path1=[DataDir filesep 'time_shift_results'];
    movefile(path1,'..');

    disp('Time shift analysis is done !');




function pre_checkbox_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The checkboxes of preprocessing  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if get(hObject,'value')
    set(handles.d2n_checkbox,'value',1);
    set(handles.RTP_checkbox,'value',1);
    set(handles.ST_checkbox,'value',1);
    set(handles.realign_checkbox,'value',1);
    set(handles.smooth_checkbox,'value',1);
    set(handles.detrend_checkbox,'value',1);
    set(handles.reg_checkbox,'value',1);
    set(handles.filter_checkbox,'value',1);
    set(handles.AM_checkbox,'value',1);
else
   
    set(handles.d2n_checkbox,'value',0);
    set(handles.RTP_checkbox,'value',0);
    set(handles.ST_checkbox,'value',0);
    set(handles.realign_checkbox,'value',0);
    set(handles.smooth_checkbox,'value',0);
    set(handles.detrend_checkbox,'value',0);
    set(handles.reg_checkbox,'value',0);
    set(handles.filter_checkbox,'value',0);
    set(handles.AM_checkbox,'value',0);
end

    

function [tc] = V_get_tc_roi(V, ROI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to read the original time courses in
% the mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mask = ROI;
Info_mask = spm_vol(Mask);
[Ymask xyz] = spm_read_vols(Info_mask);
Ymask = round(Ymask);
    ind = find(Ymask(:)>0);
    [I,J,K] = ind2sub(size(Ymask),ind);
    XYZ = [I J K]';
    XYZ(4,:) = 1;
    VY = spm_get_data(V,XYZ);
    tc= VY;
return   



function R = V_Cor_Pair(X, Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to calculate Pair-correlation between two sample matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[matX_1, matX_2] = size(X);
[matY_1, matY_2] = size(Y);

if (matX_1 ~= matY_1) || (matX_2 ~= matY_2)
    disp('Error: The two matices do not have the same size!')
else
    X = (X - repmat(mean(X, 2), 1, matX_2))./repmat(std(X, 0, 2), 1, matX_2);
    Y = (Y - repmat(mean(Y, 2), 1, matX_2))./repmat(std(Y, 0, 2), 1, matX_2);
    R = dot(X,Y,2) / (matX_2 - 1);
end



function run_btn_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RUN (preprocessing or TSA)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 handles=guidata(gcf);
 DirectoryPath=get(handles.Directory_edit,'string');
 Workdirpath=inpath_Misc(DirectoryPath,'GetParentPath');

IsPreprocessing=get(handles.pre_checkbox,'value');
IsTSA=get(handles.tsa_checkbox,'value'); 

f1=@pre_btn_Callback;
f2=@tsa_btn_Callback;

if 1==IsPreprocessing
    if 1==IsTSA
        f1(handles.pre_btn,eventdata,handles);     
        f2(handles.tsa_btn,eventdata,handles);
    else
        f1(handles.pre_btn,eventdata,handles);
    end
    
else
    if 1==IsTSA
        f2(handles.tsa_btn,eventdata,handles);
    end
end



function popupmenu3_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The popup menu of Utilities. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

select=get(hObject,'value');
switch select
    case 1

    case 2
        run automask.m
    case 3
        run coregister.m
 
end



function ref_M_btn_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The button used to select the reference maask during TSA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

handles=guidata(gcf);
    
%     DPATH=uigetdir(pwd); 
% 
%     if 0~=DPATH
%         cd(DPATH);
%         DirectoryPath=uigetfile( {'*.nii';'*.img';'*.hdr'},'File Selector');
%         workPath=[DPATH filesep DirectoryPath ];
%         set(handles.HM_edit,...
%             'string',workPath);
%      end

 [filename,pathname]= uigetfile(...
                                                       {'*.nii;*.img'});
    str_A=strcat(pathname,filename);
    if 0 ~= str_A
        cd(pathname);
        set(handles.HM_edit,...
              'string',str_A);
    end

    
    
    
    
    
    
    
    
    
    
    
    
    
    

% ----------------------------------------------------------------------------------
% OTHERS


function TR_edit_Callback(hObject, eventdata, handles)
% hObject    handle to TR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function TR_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ref_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ref_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function ref_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function slice_ord_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slice_num_edit_Callback(hObject, eventdata, handles)
% hObject    handle to slice_num_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function slice_num_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function RTP_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function low_edit_Callback(hObject, eventdata, handles)
% hObject    handle to low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function low_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function high_edit_Callback(hObject, eventdata, handles)
% hObject    handle to high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function high_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function smooth_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SR_edit_Callback(hObject, eventdata, handles)
% hObject    handle to SR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function SR_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uti_menu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in tsa_checkbox.
function tsa_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to tsa_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in d2n_checkbox.
function d2n_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to d2n_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in ST_checkbox.
function ST_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ST_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function RTP_edit_Callback(hObject, eventdata, handles)
% hObject    handle to RTP_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function smooth_edit_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function HM_edit_Callback(hObject, eventdata, handles)
% hObject    handle to HM_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function HM_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function Directory_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
