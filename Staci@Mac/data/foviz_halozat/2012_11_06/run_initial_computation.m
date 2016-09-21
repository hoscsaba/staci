clear all
clear path

global spr fname PREFIX FOVIZ_DFILE_NAME

idoszak_max=47;
is_mac=0;
spr='alapzonamodell9.spr';
FOVIZ_DFILE_NAME = 'FV_input_2012_11_06_barlang_nelkul.csv';

run_foviz_uzem(idoszak_max,is_mac);
%postprocess_opti(0,idoszak_max);