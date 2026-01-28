fileNames = dir("WS/");
fileNames = {fileNames.name};
fileNames = fileNames(4:end);
figFolder = fullfile(pwd,'Figures');
numFiles = length(fileNames);
for fileIndex = 1:numFiles
    currentName = fullfile("WS/",fileNames{fileIndex});
    load(currentName);
    numFiles = length(fileNames);
    plotmarkers.figFolder = figFolder;
    save(currentName,'plotmarkers','data','raw_data','pat_char');
end

%{
load WS/01PostVal_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/01PostVal_WS_temp.mat
clear all

load WS/01PostVal_WS.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/01PostVal_WS.mat
clear all

load WS/01PostVal2_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/01PostVal2_WS_temp.mat
clear all

load WS/01PostVal2_WS.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/01PostVal2_WS.mat
clear all

load WS/01PostVal2_god_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/01PostVal2_god_WS_temp.mat
clear all

load WS/01PreVal1_ok_kvalitet_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/01PreVal1_ok_kvalitet_WS_temp.mat
clear all

load WS/02PostVal1_OK_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/02PostVal1_OK_WS_temp.mat
clear all

load WS/02PreVal1_god_BT_mid_HR_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/02PreVal1_god_BT_mid_HR_WS_temp.mat
clear all

load WS/03PostVal2_Semidaarligefaser_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/03PostVal2_Semidaarligefaser_WS_temp.mat
clear all

load WS/03PrePracticeTydeligeFaser_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/03PrePracticeTydeligeFaser_WS_temp.mat
clear all

load WS/04PostVal2_fine_faser_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/04PostVal2_fine_faser_WS_temp.mat
clear all

load WS/04PreVal2_nogenlunde_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/04PreVal2_nogenlunde_WS_temp.mat
clear all

load WS/05PostVal1_lidt_spoejs_start_men_ellers_gode_faser_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/05PostVal1_lidt_spoejs_start_men_ellers_gode_faser_WS_temp.mat
clear all

load WS/05PreVal1_fin_kvali_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/05PreVal1_fin_kvali_WS_temp.mat
clear all

load WS/06PostVal2_paene_faser_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/06PostVal2_paene_faser_WS_temp.mat
clear all

load WS/06PreVal1_ikke_saa_tydelige_faser_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/06PreVal1_ikke_saa_tydelige_faser_WS_temp.mat
clear all

load WS/07PostVal3_rimelig_god_kvalitet_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/07PostVal3_rimelig_god_kvalitet_WS_temp.mat
clear all

load WS/07PreVal1_semi_daarlig_kval_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/07PreVal1_semi_daarlig_kval_WS_temp.mat
clear all

load WS/08PostVal2_god_kval_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/08PostVal2_god_kval_WS_temp.mat
clear all

load WS/08PreVal1_semi_daarlig_kval_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/08PreVal1_semi_daarlig_kval_WS_temp.mat
clear all

load WS/09PostVal_Val1_ok_kval_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/09PostVal_Val1_ok_kval_WS_temp.mat
clear all

load WS/09PreVal2_god_kvalitet_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/09PreVal2_god_kvalitet_WS_temp.mat
clear all

load WS/10PostVal1_OK_men_BT_overstiger_200_og_forsvinder_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/10PostVal1_OK_men_BT_overstiger_200_og_forsvinder_WS_temp.mat
clear all

load WS/10PreVal2_ok_kval_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/10PreVal2_ok_kval_WS_temp.mat
clear all

load WS/14PostVal1_Meget_god_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/14PostVal1_Meget_god_WS_temp.mat
clear all

load WS/Subject_1_WS_temp.mat
figFolder = fullfile(pwd,'Figures');
plotmarkers.figFolder = figFolder;
save WS/Subject_1_WS_temp.mat
clear all
%}



