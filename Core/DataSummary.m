function [T_xl] = DataSummary(list)
% function DataSummary
% Input: list of quantities to save to data file (could be mult. patients)
% Output: excel spreadsheet with computed quantities
% Description: mat2excel Writes markers from mat files to an excel spreadsheet

cd .. 
% Preallocating info arrays
Info_xl = strings(length(list),5);
Info_csv = zeros(length(list),5);
ModelParameters = zeros(length(list),10);

for i = 1:length(list)
    patient   = list{i};
    patient   = patient(1:end-4);

    % Creating string info array for xl
    load(strcat('WS/',patient,'_WS.mat'))

    PatientID   = pat_char.PatientID;
    Age         = pat_char.age;
    Sex         = pat_char.sex;
    Height      = pat_char.height;
    Weight      = pat_char.weight;

    Info_xl(i,1)    = PatientID;
    Info_xl(i,2)    = Age; 
    Info_xl(i,3)    = Sex;
    Info_xl(i,4)    = Height;
    Info_xl(i,5)    = Weight;

    % Creating number info array for csv
    Info_csv(i,1)    = PatientID;
    Info_csv(i,2)    = Age;        %%% add patient ID number to popup window

    % Changing sex to number value for the csv(m:1, f:2)                                                                 for loop section
    if Info_xl(i,3) == 'm'
        Info_csv(i,3) = 1;
    else 
        Info_csv(i,3) = 2;
    end

    Info_csv(i,4)    = Height;
    Info_csv(i,5)    = Weight;

    % make different info table for csv (only num) (T_xl, T_csv)
    
    T_ClinicalRatios(i,:) = readtable(strcat('Markers/',patient,'_markers.xlsx'));
    
    load(strcat('Optimized/',patient,'_opt.mat'),'optpars')
    optpars = exp(optpars);
    ModelParameters(i,1)  = optpars(19);  %xi_w  patient specific
    ModelParameters(i,2)  = optpars(20);  %xi_p  patient specific
    ModelParameters(i,3)  = optpars(21);  %xi_r  patient specific
    ModelParameters(i,4)  = optpars(22);  %xi_s  patient specific
    ModelParameters(i,5)  = optpars(23);  %H_I   patient specific
    ModelParameters(i,6)  = optpars(24);  %H_p   estimated
    ModelParameters(i,7)  = optpars(25);  %H_r   estimated
    ModelParameters(i,8)  = optpars(26);  %H_s   estimated
    ModelParameters(i,9)  = optpars(27);  %H_Ia  patient specific
    ModelParameters(i,10) = optpars(28);  %tchar estimated
end

% T_Info_csv might be redundant bc it has headers which are ultimately
%    ignored but keeping for now so that T_csv works
T_ClinicalRatios = removevars(T_ClinicalRatios,{'Var1'});
T_Info_xl = array2table(Info_xl, 'VariableNames', {'Patient', 'Age', 'Sex','Height','Weight'});
T_Info_csv = array2table(Info_csv, 'VariableNames', {'Patient', 'Age', 'Sex','Height','Weight'});
T_ModelParameters = array2table(ModelParameters, 'VariableNames', ...
    {'xi_w', 'xi_p','xi_r','xi_s','H_I','H_p','H_r','H_s','HIa','tchar'});
T_xl = [T_Info_xl, T_ClinicalRatios, T_ModelParameters];
T_csv = [T_Info_csv, T_ClinicalRatios, T_ModelParameters];

% Input window code
prompt = {'Enter data summary file name (ex: FileName)'};
dlgtitle = 'Save Data';
dims = [1 36];
definput = {'FileName','hsv'};
answer = inputdlg(prompt,dlgtitle,dims,definput);

% Titles for excel and csv files
xl_answer = append(answer{1},'.xlsx');
csv_answer = append(answer{1},'.csv');

if isempty(answer) == 1
    return
end

% Saves file to Excel
if isfile(xl_answer)
    disp('file exist')
    quest = {'Add or overwrite existing data?'};
    dlgtitle = 'Save Data';
    answer2 = questdlg(quest,dlgtitle,'Add','Overwrite','Add');
    if strcmp(answer2, 'Add') == 1
        writetable(T_xl,xl_answer,'WriteMode','Append','AutoFitWidth',false)
        writetable(T_csv,csv_answer,'WriteMode','Append') 
    else
        writetable(T_xl,xl_answer, 'WriteMode','Overwrite')
        writetable(T_csv,csv_answer, 'WriteMode','Overwrite')
    end
else
    disp('file does not exist')
    writetable(T_xl,xl_answer, 'WriteMode','Overwrite')
    writetable(T_csv,csv_answer, 'WriteMode','Overwrite')
end

cd Core

end
