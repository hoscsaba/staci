function out = sort_data(data)

% ADATOK SZORTIROZASA
% 1 - Staci , FV folyamatvaltozo, Adattipus
%  2 - #med_barlang, K_KB_GELL_BARL_01KAM_HS_M-AVGHALF, Vizszint (m)
%  3 - #med_barlang, -, Zarallapot
%  4 - #med_cinkota, K_EP_CINK_MED_01KAM_HS_M-AVGHALF, Vizszint (m)
%  5 - #med_cinkota, -, Zarallapot
%  6 - #med_gilice, K_DP_GILI_UMED_KAM_HS_M-AVGHALF, Vizszint (m)
%  7 - #med_gilice, -, Zarallapot
%  8 - #med_kob_regi, K_DP_KOBA_RMED_01KAM_HS_M-AVGHALF, Vizszint (m)
%  9 - #med_kob_regi, -, Zarallapot
%  10 - #med_kob_uj, K_DP_KOBA_UMED_02KAM_HS_M-AVGHALF, Vizszint (m)
%  11 - #med_kob_uj, -, Zarallapot
%  12 - #med_krisztina, K_KB_KRIS_UMED_01KAM_HS_M-AVGHALF, Vizszint (m)
%  13 - #med_krisztina, -, Zarallapot
%  14 - #med_rakossz, K_EP_RSZM_MED_KAM_HS_M-AVGHALF, Vizszint (m)
%  15 - #med_rakossz, -, Zarallapot
%  16 - #med_sanc, K_KB_SANC_MED_02KAM_HS_M-AVGHALF, Vizszint (m)
%  17 - #med_sanc, -, Zarallapot
%  18 - Alapzona, K_K_KOZP_PEST_BUDA_QOSSZ_SM-AVGHALF-K_DB_BORS_GH_Q_M-AVGHALF-Krisztinak-KEPE, Fogyasztas (m3/h)
%  19 - KeletPest, K_K_KOZP_KEPE_FELSOZONA_QOSSZ_SM-AVGHALF, Fogyasztas (m3/h)
%  20 - A_FELSOJH_GH, K_EB_BUJL_TERM_Q_SM-AVGHALF, Betap Q (m3/h)
%  21 - A_Z20-RADNOTI, K_EB_MARG_RGH_QTERM_SM-AVGHALF, Betap Q (m3/h)
%  22 - A_BEKAS_PUM, E_BM_VT_BM_ROKA_QOS_SM-AVGHALF, Betap Q (m3/h)
%  23 - A_KM4_Q1600, E_BP_VT_NNYGH_Q_SM-AVGHALF, Betap Q (m3/h)
%  24 - A_KPM1-3, E_KM_VT_12GH_Q_SM-AVGHALF, Betap Q (m3/h)
%  25 - A_CSEPEL_BF, D_CS_VT_GH_Q_SM-AVGHALF, Betap Q (m3/h)
%  26 - A_BUDAORSI_UT, K_DB_BORS_GH_Q_M-AVGHALF, Tovabbemelt Q (m3/h)
%  27 - A_KRISZTINA_UJ, K_KB_KRIS_DGH_Q_SM-AVGHALF +K_KB_KRIS_SASLIP_QOSSZ_SM-AVGHALF+K_KB_KRIS_VGH_Q_SM-AVGHALF, Tovabbemelt Q (m3/h)
%  28 - K_KOBANYA_GH, K_DP_KOBA_PILL_Q_SM-AVGHALF, Tovabbemelt Q (m3/h)
%  29 - K_RSZM_GH, K_EP_RSZM_GH_Q_M-AVGHALF, Tovabbemelt Q (m3/h)
%  30 - K_GILICE_UJ_KEPE_GH, K_DP_GILI_KEPE_Q_SM-AVGHALF, Tovabbemelt Q (m3/h)

%  SPECI CSOMOPONTOK
idx_spec_nodes=1;

for j=1:length(data)
	%clc
	j
	data{j}
	fprintf('\n %d - %s, %s, %s',j,data{j}{1},data{j}{2},data{j}{3});
end

% MEDENCEK

for j=1:7
	out.med{j}.name=data{2*j}{1};
	for i=1:48
		out.med{j}.level(i)=str2num(data{2*j}{3+i});
		if strcmp(data{2*j+1}{3+i},'kizarva')
			out.med{j}.status(i)=0;
		elseif strcmp(data{2*j+1}{3+i},'toltodik')
			out.med{j}.status(i)=1;
		elseif strcmp(data{2*j+1}{3+i},'urul') || strcmp(data{2*j+1}{3+i},'orol')
			out.med{j}.status(i)=-1;
		elseif strcmp(data{2*j+1}{3+i},'Ellennyomo')
			out.med{j}.status(i)=2;
		else
			data{2*j+1}{3+i}
			error('Ismeretlen medence allapot!!!');
		end
	end
end

% ALAPZONA FOGYASZTAS
for i=1:48
	out.alapzona(i) =str2num(data{18-2}{3+i});
	out.keletpest(i)=str2num(data{19-2}{3+i});
end

% BETAPOK
for j=1:6
	out.betap{j}.name=data{19+j-2}{1};
	out.spec_nodes{idx_spec_nodes}=data{19+j-2}{1};
	idx_spec_nodes=idx_spec_nodes+1;
	for i=1:48
		out.betap{j}.Q(i) =str2num(data{19+j-2}{3+i});
	end
end

% TOVABBEMELESEK
for j=1:8
	out.tovabbemeles{j}.name=data{25+j-2}{1};
    out.tovabbemeles{j}.is_internal=0;

	% Tovabbemelo GH eseten masik zona csatlakozasi pontja
	out.tovabbemeles{j}.twin='';
	out.tovabbemeles{j}.mul=1;
	if strcmp(out.tovabbemeles{j}.name,'K_KOBANYA_GH')==1
		out.tovabbemeles{j}.twin='NODE317';
		out.tovabbemeles{j}.mul=-1;
	end

	if strcmp(out.tovabbemeles{j}.name,'K_RSZM_GH')==1
		out.tovabbemeles{j}.twin='NODE293';
		out.tovabbemeles{j}.mul=-1;
	end

	if strcmp(out.tovabbemeles{j}.name,'K_GILICE_UJ_KEPE_GH')==1
		out.tovabbemeles{j}.twin='NODE305';
		out.tovabbemeles{j}.mul=-1;
    end
    
    if strcmp(out.tovabbemeles{j}.name,'NODE320')==1
		out.tovabbemeles{j}.twin='NODE688';
		out.tovabbemeles{j}.mul=1;
        out.tovabbemeles{j}.is_internal=1;
	end

	out.spec_nodes{idx_spec_nodes}=data{25+j-2}{1};
	idx_spec_nodes=idx_spec_nodes+1;
	for i=1:48
		out.tovabbemeles{j}.Q(i) =str2num(data{25+j-2}{3+i});
	end
end

% OSSZESITES
for i=1:48
	out.sumQki(i)=out.alapzona(i)+out.keletpest(i);
	for j=1:length(out.tovabbemeles)
		out.sumQki(i)=out.sumQki(i)+out.tovabbemeles{j}.Q(i);
	end

	out.sumQbe(i)=0;
	for j=1:length(out.betap)
		out.sumQbe(i)=out.sumQbe(i)+out.betap{j}.Q(i);
	end
end

% MEDENCE ZAR ADATOK
	% out.med.status
	% kizarva  = 0
	% toltodik = 1
	% urul = -1
	% ellennyomo = 2


% valve_status_vals:
% sor: status
% oszlop: valve

%  2 - #med_barlang
% out.med{1}.staci_name = 'POOL584';
% out.med{1}.A = 5000;
% out.med{1}.valves = {'VALVE572'};
% out.med{1}.valve_status_name = {'kizarva'};
% out.med{1}.valve_status_idx  = [0];
% out.med{1}.valve_status_vals = [0];

%  4 - #med_cinkota
med_idx=1;
out.med{med_idx}.staci_name = 'POOL288';
out.med{med_idx}.A = 1670;
out.med{med_idx}.valves = {};
out.med{med_idx}.valve_status_name = [];
out.med{med_idx}.valve_status_idx  = [];
out.med{med_idx}.valve_status_vals = [];

%  6 - #med_gilice
med_idx=med_idx+1;
out.med{med_idx}.staci_name = 'POOL315';
out.med{med_idx}.A = 2783;
out.med{med_idx}.valves = {'VALVE312','VALVE306','VALVE309'};
out.med{med_idx}.valve_status_name = {'toltodik','urul','kizarva'};
out.med{med_idx}.valve_status_idx  = [1 -1 0];
out.med{med_idx}.valve_status_vals = [0 1 1
1 0 0
0 1 0];

%  8 - #med_kob_regi
med_idx=med_idx+1;
out.med{med_idx}.staci_name = 'POOL327';
out.med{med_idx}.A = 2750/2; % Egyik kamra ki van zárva!
out.med{med_idx}.valves = {'VALVE321'};
out.med{med_idx}.valve_status_name = {'kizarva','urul','toltodik'};
out.med{med_idx}.valve_status_idx  = [0 -1 1];
out.med{med_idx}.valve_status_vals = [0
0
1];

%  10 - #med_kob_uj
med_idx=med_idx+1;
out.med{med_idx}.staci_name = 'POOL324';
out.med{med_idx}.A = 3333;
out.med{med_idx}.valves = {'VALVE617','VALVE318'};
out.med{med_idx}.valve_status_name = {'toltodik','urul','kizarva'}; %volt {'kizarva','urul','toltodik'}; by hp 03.25
out.med{med_idx}.valve_status_idx  = [1 -1 0];
out.med{med_idx}.valve_status_vals = [1 1
1 0
0 1];

%  12 - #med_krisztina
med_idx=med_idx+1;
out.med{med_idx}.staci_name = 'POOL276';
out.med{med_idx}.A = 3767;
out.med{med_idx}.valves = {'VALVE273','VALVE270','VALVE267'};
out.med{med_idx}.valve_status_name = {'toltodik','urul','kizarva'};
out.med{med_idx}.valve_status_idx  = [1 -1 0];
out.med{med_idx}.valve_status_vals = [0 1 1 %1 1 0 %Mod by Haraszti Péter 06.01.27.
1 0 0 % 1 0 0 Mod by HP+HCs 06.02.04.
0 0 1];

%  14 - #med_rakossz
med_idx=med_idx+1;
out.med{med_idx}.staci_name = 'POOL303';
out.med{med_idx}.A = 1667;
out.med{med_idx}.valves = {'VALVE297','VALVE300','VALVE294'};
out.med{med_idx}.valve_status_name = {'toltodik','urul','kizarva'};
out.med{med_idx}.valve_status_idx  = [1 -1 0];
out.med{med_idx}.valve_status_vals = [1 0 1
0 1 0
0 0 1];

%  16 - #med_sanc
med_idx=med_idx+1;
out.med{med_idx}.staci_name = 'POOL336';
out.med{med_idx}.A = 10000;
out.med{med_idx}.valves = {'VALVE686','VALVE588','VALVE575'};
out.med{med_idx}.valve_status_name = {'ellennyomo','kizarva'};
out.med{med_idx}.valve_status_idx  = [2 0];
out.med{med_idx}.valve_status_vals = [1 1 1
    0 1 1];
end


function plot_input_data(out)
% ABRAK
figure(1)
time=(0:1:length(out.med{1}.level)-1)/2;
for j=1:8
	subplot(4,2,j)
	plot(time,out.med{j}.level), hold on
	bar(time,out.med{j}.status,'r'), hold off
	grid on
	xlabel('t [h]'), ylabel('h [m]')
	title(out.med{j}.name,'Interpreter','none')
	xlim([0 24])
end

figure(2)
subplot(2,1,1)
plot(time,out.alapzona/1000,'r',time,out.keletpest/1000,'b')
legend('Alapzona','Kelet-Pest')
xlabel('t [h]'), ylabel('Q, [ezer m^3/h]')
ylim([0 5])

subplot(2,1,2)
plot(time,out.sumQbe/1000,'r',time,out.sumQki/1000,'b')
legend('Betap','Fogyasztas+tovabbemeles')
xlabel('t [h]'), ylabel('Q, [ezer m^3/h]')
% ylim([0 5000])

end
