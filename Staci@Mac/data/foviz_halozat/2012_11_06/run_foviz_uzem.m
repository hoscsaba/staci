function run_foviz_uzem(idoszak_max,is_mac)
	
	clc

	global logfname spr Q_A Q_K DEBUG
	global fname spr FOVIZ_DFILE_NAME

	path1 = getenv('PATH')
	if is_mac==1
		path1 = [path1 ':/Users/hoscsaba/Dropbox/staci/kod/new_staci/Staci@Mac'];
	else
		path1 = [path1 ':/Users/hoscsaba/Dropbox/staci/kod/new_staci/Staci@Wylie'];
	end
	setenv('PATH', path1);

	DEBUG=1;
	spr='alapzonamodell9.spr';
	fname = FOVIZ_DFILE_NAME;
	logfname = 'log_run_foviz_uzem.txt';
	fp=fopen(logfname,'w');
	fclose(fp);

	data = read_dfile(fname);
	data = sort_data(data);

	Q_A=get_zone_consumptions(spr,data,'A_');
	Q_K=get_zone_consumptions(spr,data,'K_');

	fp=fopen(logfname,'a');
	mfprintf([1 fp],'\n\n   BUILDING SOLVER INPUT FILES...\n======================================');
	fclose(fp);

	for idoszak=1:idoszak_max
		fp=fopen(logfname,'a');
		mfprintf([1 fp],'\n\n   IDOSZAK: %g\n-------------------------------',idoszak);
		fclose(fp);

		% Filenev beallitasa
		spr_new=[fname(1:end-4),'_',num2str(idoszak),'.spr'];
		copyfile(spr,spr_new);

		% Fogyasztasok beallitasa
		set_all_consumptions(idoszak,data,spr_new);

		check_zone_balance(idoszak,spr_new,data,'A_');
		check_zone_balance(idoszak,spr_new,data,'K_');
		
		% Medence zarallapotok
		set_pools(idoszak,data,spr_new);
	end
end

%------------------------------------------------
function [H_HDS,H_FV,nodenames]=get_spec_pres(idoszak,data,spr_new)

	global logfname spr

	fp=fopen(logfname,'a');
	mfprintf([1 fp],'\n\n   Nyomasok ellenorzese:');

	data=read_dfile('FV_input_2012_11_06_eredmenyek_foviz.csv');
	save('tmp.mat','data');
	kk=1;
	for i=2:length(data)
		FV_ID=data{i}{1};
		Staci_ID=data{i}{2};
		if length(Staci_ID)>2
			H_FV(kk) = str2num(data{i}{4+idoszak});
			H_HDS(kk) = get_node_head(spr_new,Staci_ID);
			nodenames{kk}=Staci_ID;
			mfprintf([1 fp],'\n\t FV: %s, Staci alias: %s, nyomas: FV : %g vom, HDS : %g vom',FV_ID,Staci_ID,H_FV(kk),H_HDS(kk));
			kk=kk+1;
		end
	end

end

%------------------------------------------------
function [QpoolHDS,QpoolFVM]=get_pool_flow_rates(idoszak,data,fname)

	global logfname spr

	fp=fopen(logfname,'a');
	mfprintf([1 fp],'\n\n   Medencek staci terfogatarama:');

	for i=1:length(data.med)
		mfprintf([1 fp],'\n     %s :',add_whitespace(data.med{i}.name,14));

		QpoolHDS(i)=get_edge_flowrate(fname,data.med{i}.staci_name);
		QpoolFVM(i)=data.med{i}.A*(data.med{i}.level(idoszak+1)-data.med{i}.level(idoszak))*2;
		mfprintf([1 fp],'  %+6.2f m3/h, Foviz adat: %+6.2f',QpoolHDS(i),QpoolFVM(i));
		mfprintf([1 fp],'\t  %+6.2f = ( %g - %g ) * %g * 2',QpoolFVM(i),...
			data.med{i}.level(idoszak+1),data.med{i}.level(idoszak),data.med{i}.A);

	end
	fclose(fp);
end

%------------------------------------------------
function set_pools(idoszak,data,fname)

	global logfname spr

	fp=fopen(logfname,'a');
	mfprintf([1 fp],'\n\n   Medencek beallitasa:');

	% out.med.status
	% kizarva  = 0
	% toltodik = 1
	% urul = -1
	% ellennyomo = 2
	
	for i=1:length(data.med)
		mfprintf([1 fp],'\n   %s :',data.med{i}.name);

		set_pool_water_level(fname,data.med{i}.staci_name,data.med{i}.level(idoszak));
		mfprintf([1 fp],' vizszint = %g m, ',data.med{i}.level(idoszak));

		if ~isempty(data.med{i}.valve_status_idx)			
			status_idx = find(data.med{i}.valve_status_idx==data.med{i}.status(idoszak));
			if isempty(status_idx)
				error('WTF???');
			end
			mfprintf([1 fp],' %s',data.med{i}.valve_status_name{status_idx});

			for j=1:length(data.med{i}.valves)
				opening = abs(data.med{i}.valve_status_vals(status_idx,j)-1.)*100.;
				if data.med{i}.valve_status_vals(status_idx,j)==1
					str='nyitva';
				elseif data.med{i}.valve_status_vals(status_idx,j)==0
					str='zarva ';
				else
					error('VALVE WTF????')
				end
				mfprintf([1 fp],'\n\t\t %s : %s, allas = %g %%',data.med{i}.valves{j},str,opening);
				set_valve_opening(fname,data.med{i}.valves{j},opening);		
				%set_valve_opening(fname,data.med{i}.valves{j},0);			
			end
		end
	end
	fclose(fp);
end

%------------------------------------------------
function is_ok = run_staci(fname)
	global logfname spr

	fp=fopen(logfname);
	fname_rrs=[fname,'.rrs'];
	if exist(fname_rrs,'file')==2
		mfprintf([1 fp],'\n\n Deleting %s ...',fname_rrs);
		delete(fname_rrs);
	end

	staci_cmd = sprintf('system(''new_staci -s %s'')',fname);
	mfprintf([1 fp],'\n\n Running %s ...',staci_cmd);
	
	fclose(fp);

	res=strtrim(evalc(staci_cmd));

	fp=fopen(logfname);
	mfprintf([1 fp],'%s',res);
	fclose(fp);

	is_ok=strcmp(strtrim(fileread(fname_rrs)),'ok');
end

%------------------------------------------------
function set_pool_water_level(fname,poolname,value)
	staci_cmd = sprintf('system(''new_staci -m %s -e %s -p water_level -n %g -o %s'')',...
		fname,poolname,value,'tmp.spr');
	res=strtrim(evalc(staci_cmd));

	copyfile('tmp.spr',fname);
	delete('tmp.spr');
end

%------------------------------------------------
function set_valve_opening(fname,valvename,value)
	staci_cmd = sprintf('system(''new_staci -m %s -e %s -p position -n %g -o %s'')',...
		fname,valvename,value,'tmp.spr');
	res=strtrim(evalc(staci_cmd));

	copyfile('tmp.spr',fname);
	delete('tmp.spr');
end

%------------------------------------------------
function set_all_consumptions(idoszak,data,fname)

	global logfname spr Q_A Q_K

	% Zonafogyasztasok szamitasa es beallitasa
	Q_A_kivant = data.alapzona(idoszak);
	Q_K_kivant = data.keletpest(idoszak);

	fp=fopen(logfname,'a');
	mfprintf([1 fp],'\n   kivant fogyasztasok: Alapzona   : %g m3/h',Q_A_kivant);
	mfprintf([1 fp],'\n                        Kelet Pest : %g m3/h',Q_K_kivant);
	fclose(fp);
Q_A_kivant
Q_A
	set_zone_consumptions(fname,data,'A_',Q_A_kivant/Q_A);
	set_zone_consumptions(fname,data,'K_',Q_K_kivant/Q_K);

	% ELLENORZES
	%Q_A_uj=get_zone_consumptions(fname,data,'A_');
	%Q_K_uj=get_zone_consumptions(fname,data,'K_');

	% SPECIALIS CSOMOPONTOK
	fp=fopen(logfname,'a');
	mfprintf([1 fp],'\n   Specialis csompontok:');

	node_list=importdata('node_list.csv');

% BETAPOK
for j=1:length(data.betap)
	megvan=0;
	for i=1:length(node_list)
		% node_list{i} 
		% data.betap{j}.name
		% fprintf('\n====================');
		% pause
		if strcmp(node_list{i},data.betap{j}.name)==1
			megvan=1;
			Qcsp=-data.betap{j}.Q(idoszak);
			set_node_consumption(fname,node_list{i},Qcsp);

			node_name = add_whitespace(node_list{i},22);
			mfprintf([1 fp],'\n \t %s:\t %6.2f m3/h (betap)',node_name,Qcsp);
		end
	end
	if (megvan==0)
		fp=fopen(logfname,'a');
		mfprintf([1 fp],'\n   Nem talalom ezt a specialis csomopont: %s\n',data.betap{j}.name);
		fclose(fp);
		error('PANIK!');
	end
end

% TOVABBEMELESEK
for j=1:length(data.tovabbemeles)
	megvan=0;
	for i=1:length(node_list)
		if strcmp(node_list{i},data.tovabbemeles{j}.name)==1
			megvan=1;
			Qcsp=data.tovabbemeles{j}.Q(idoszak)*data.tovabbemeles{j}.mul;
			set_node_consumption(fname,node_list{i},Qcsp);

			node_name = add_whitespace(node_list{i},22);
			mfprintf([1 fp],'\n \t %s:\t %6.2f m3/h (tovabbemeles)',node_name,Qcsp);

			if ~isempty(data.tovabbemeles{j}.twin)
				set_node_consumption(fname,data.tovabbemeles{j}.twin,-Qcsp);
				mfprintf([1 fp],' + ikercsp: %s -> %6.2f m3/h',data.tovabbemeles{j}.twin,-Qcsp);
			end
		end
	end
	if (megvan==0)
		fp=fopen(logfname,'a');
		mfprintf([1 fp],'\n   Nem talalom ezt a specialis csomopont: %s\n',data.betap{j}.name);
		fclose(fp);
		error('PANIK!');
	end
end

fclose(fp);


end

%------------------------------------------------
function check_zone_balance(idoszak,fname,data,prefix)

	global logfname 

	fp=fopen(logfname,'a');
	mfprintf([1 fp],'\n\n Vizmerleg ellenorzese, adatfajl: %s\n',fname);

	Qfogy=get_zone_consumptions(fname,data,prefix);
	Qzona=-Qfogy;

	if strcmp(prefix,'A_')==1
		mfprintf([1 fp],'\n\t Alapzona fogyasztas    : \t%g m3/h\n',Qzona);

	% BETAPOK
	for i=1:length(data.betap)
		Qin=-get_node_consumption(fname,data.betap{i}.name);
		Qzona=Qzona+Qin;
		node_name = add_whitespace(data.betap{i}.name,22);
		mfprintf([1 fp],'\n \t %s:\t %6.1f m3/h (betap), Qzona=%6.1f m3/h',node_name,Qin,Qzona);
	end

	% TOVABBEMELESEK
	for i=1:length(data.tovabbemeles)
		if data.tovabbemeles{i}.is_internal==0
			Qout=get_node_consumption(fname,data.tovabbemeles{i}.name)*data.tovabbemeles{i}.mul;
			Qzona=Qzona-Qout;
			node_name = add_whitespace(data.tovabbemeles{i}.name,22);
			mfprintf([1 fp],'\n \t %s:\t %6.1f m3/h (tovabbemeles), Qzona=%6.1f m3/h',node_name,Qout,Qzona);
		else
			node_name = add_whitespace(data.tovabbemeles{i}.name,22);
			Qout=get_node_consumption(fname,data.tovabbemeles{i}.name)*data.tovabbemeles{i}.mul;
			mfprintf([1 fp],'\n \t %s:\t %6.1f m3/h (belso tovabbemeles, figyelmen kivul hagyva)',node_name,Qout);
		end
	end

	% MEDENCEK
	mfprintf([1 fp],'\n');
	for i=1:length(data.med)
		mfprintf([1 fp],'\n \t %s:\t',add_whitespace(data.med{i}.name,22));
		if (strcmp('#med_cinkota',data.med{i}.name)==0)
			Qpool=data.med{i}.A*(data.med{i}.level(idoszak+1)-data.med{i}.level(idoszak))*2;
			Qzona=Qzona-Qpool;
			mfprintf([1 fp],'  %+6.1f m3/h, Qzona: %+6.1f',Qpool,Qzona);
			mfprintf([1 fp],'\t  %+6.2f = ( %g - %g ) * %g * 2',Qpool,...
				data.med{i}.level(idoszak+1),data.med{i}.level(idoszak),data.med{i}.A);
		else
			mfprintf([1 fp],'  atugrom, nem Kelet-Pestzona.');
		end
	end

elseif strcmp(prefix,'K_')==1
	mfprintf([1 fp],'\n\t KELET-PEST fogyasztas  : \t%g m3/h\n',Qzona);

	% TOVABBEMELESEK = BETERMELES
	for i=1:length(data.tovabbemeles)
		if strncmp(data.tovabbemeles{i}.name,'K_',2)==1 
			Qout=get_node_consumption(fname,data.tovabbemeles{i}.twin);
			Qzona=Qzona+Qout;
			node_name = add_whitespace(data.tovabbemeles{i}.twin,22);
			mfprintf([1 fp],'\n \t %s:\t %6.1f m3/h (betap Kelet-Pestre), Qzona=%6.1f m3/h',node_name,Qout,Qzona);
		end
	end

	% MEDENCEK
	mfprintf([1 fp],'\n');
	for i=1:length(data.med)
		mfprintf([1 fp],'\n \t %s:\t',add_whitespace(data.med{i}.name,22));
		if (strcmp('#med_cinkota',data.med{i}.name)==1)
			Qpool=data.med{i}.A*(data.med{i}.level(idoszak+1)-data.med{i}.level(idoszak))*2;
			Qzona=Qzona-Qpool;
			mfprintf([1 fp],'  %+6.1f m3/h, Qzona: %+6.1f',Qpool,Qzona);
			mfprintf([1 fp],'\t  %+6.2f = ( %g - %g ) * %g * 2',Qpool,...
				data.med{i}.level(idoszak+1),data.med{i}.level(idoszak),data.med{i}.A);
		else
			mfprintf([1 fp],'  atugrom, nem alapzona.');
		end
	end
else
	error('Ismeretelen ZONA!!!')
end	

mfprintf([1 fp],'\n\n\t Zona vizmerleg hiba: %g m3/h (%g%%)\n\n',Qzona,abs(100*Qzona/Qfogy));


fclose(fp);

end

%------------------------------------------------
function set_zone_consumptions(fname,data,prefix,mul)

	global logfname DEBUG

	fp=fopen(logfname,'a');
	mfprintf([1 fp],'\n\n Zonafogyasztasok modositasa, adatfajl: %s, ',fname);

	if strcmp(prefix,'A_')
		mfprintf([1 fp],' Alapzona\n');
	elseif strcmp(prefix,'K_')
		mfprintf([1 fp],' Kelet-Pest\n');
	else
		error('WTF???');
	end

	node_list=importdata('node_list.csv');

	idx=strncmp(node_list,prefix,2);

	for i=1:length(idx)
		if idx(i)==1 
			Qcsp=get_node_consumption(fname,node_list{i});
			node_name = add_whitespace(node_list{i},22);
			if DEBUG>0
				mfprintf([1 fp],'\n \t %s:\t %6.2f m3/h',node_name,Qcsp);
			end
			if sum(strcmp(node_list{i},data.spec_nodes))==0
				Q=Qcsp*mul;
				set_node_consumption(fname,node_list{i},Q);
				if DEBUG>0
					mfprintf([1 fp],'\t -> %6.2f m3/h',Q);
					Q1=get_node_consumption(fname,node_name);
					mfprintf([1 fp],' (ellenorzes: Q=%g m3/h)',Q1);
				end
			else
				if DEBUG>0
					mfprintf([1 fp],'\t specialis csomopont');
				end
			end			
		end
	end
	
	fclose(fp);

end

%------------------------------------------------
function out=add_whitespace(str,max_length)
	out=str;
	if length(str<max_length)
		for kk=length(str):max_length
			out=[out,' '];
		end
	end
end

%------------------------------------------------
function set_node_consumption(fname,nodename,value)
	global logfname DEBUG
	staci_cmd = sprintf('system(''new_staci -m %s -e %s -p demand -n %g -o %s'')',...
		fname,nodename,value,'tmp.spr');
	res=strtrim(evalc(staci_cmd));

	% if DEBUG>0
	% 	fp=fopen(logfname,'a');
	% 	mfprintf([1 fp],'\n \t %s demand set to \t %6.2f m3/h',nodename,value);
	% 	fclose(fp);
	% end

	copyfile('tmp.spr',fname);
	delete('tmp.spr');
end

%------------------------------------------------
function Qcsp=get_node_consumption(fname,nodename)
	staci_cmd = sprintf('system(''new_staci -g %s -e %s -p demand'')',fname,nodename);

	res=strtrim(evalc(staci_cmd));

	idx1=strfind(res,'value');
	idx2=strfind(res,'Timing');
	[tok,rem]=strtok(res(idx1+6:idx2-1),':');
	Qcsp=str2num(strtrim(rem(2:end)));

end

%------------------------------------------------
function pcsp=get_node_head(fname,nodename)
	staci_cmd = sprintf('system(''new_staci -g %s -e %s -p head'')',fname,nodename);

	res=strtrim(evalc(staci_cmd));

	idx1=strfind(res,'value');
	idx2=strfind(res,'Timing');
	[tok,rem]=strtok(res(idx1+6:idx2-1),':');
	pcsp=str2num(strtrim(rem(2:end)));
end


%------------------------------------------------
function Qag=get_edge_flowrate(fname,edgename)
	staci_cmd = sprintf('system(''new_staci -g %s -e %s -p mass_flow_rate'')',fname,edgename);

	res=strtrim(evalc(staci_cmd));

	idx1=strfind(res,'value');
	idx2=strfind(res,'Timing');
	[tok,rem]=strtok(res(idx1+6:idx2-1),':');
	Qag=str2num(strtrim(rem(2:end)))*3.6;
end

%------------------------------------------------
function Q=get_zone_consumptions(fname,data,prefix)

	global logfname DEBUG

	fp=fopen(logfname,'a');
	mfprintf([1 fp],'\n\n Zonafogyasztasok kiszamitasa, adatfajl: %s\n',fname);

	node_list=importdata('node_list.csv');

	idx=strncmp(node_list,prefix,2);

	Q=0;

	for i=1:length(idx)
		if idx(i)==1 
			Qcsp=get_node_consumption(fname,node_list{i});
			node_name = add_whitespace(node_list{i},22);
			if DEBUG>0
				mfprintf([1 fp],'\n \t %s:\t %6.2f m3/h',node_name,Qcsp);
				%pause
			end
			if sum(strcmp(node_list{i},data.spec_nodes))==0
				Q=Q+Qcsp;
				if DEBUG>0
					mfprintf([1 fp],'\t zona : %6.2f m3/h',Q);
				end
			else
				if DEBUG>0
					mfprintf([1 fp],'\t specialis csomopont');
				end
			end			
		end
	end
	if strcmp(prefix,'A_')
		tmp='Alapzona ';
	elseif strcmp(prefix,'K_')
		tmp='Kelet-Pest ';
	else
		error('WTF???');
	end
	mfprintf([1 fp],'\n\n %s osszesen: %g m3/h\n',tmp,Q);
	fclose(fp);

end

%------------------------------------------------
function set_consumptions(data)

	node_list=importdata('node_list.csv');

	e_name=node_list(1)

	staci_cmd = sprintf('system(''new_staci -m alapmodell.spr'')');
	res=evalc(staci_cmd);

end

%------------------------------------------------
function out = read_dfile(fname)

fprintf('%s adatfajl olvasasa...',fname);

% Replace , by .

fin = fopen(fname);
fout = fopen('tmp.txt','w');

while ~feof(fin)
   s = fgetl(fin);
   s = strrep(s, ',', '.');
   fprintf(fout,'%s\n',s);
  % disp(s)
end

fclose(fin);
fclose(fout);

copyfile(fname,[fname,'.old']);
copyfile('tmp.txt',fname);

fp=fopen(fname);
line_num=1;
while ~feof(fp)
    tmp=fgetl(fp);
    i=1;
    while ~isempty(tmp)
        [out{line_num}{i},tmp]=strtok(tmp,';');
        i=i+1;
    end
    line_num=line_num+1;
end

fprintf(' kesz.\n');
fclose(fp);
end

%------------------------------------------------
function out = read_dfile2(fname)

	fprintf('%s adatfajl olvasasa...',fname);

	fp=fopen(fname);
	line_num=1;
	while ~feof(fp)
		tmp=fgetl(fp);
		i=1;
		while ~isempty(tmp)
			[out{line_num}{i},tmp]=strtok(tmp,';');
			i=i+1;
		end
		line_num=line_num+1;
	end

	fprintf(' kesz.\n');

end

%------------------------------------------------

function mfprintf(fid, varargin)
	arrayfun(@(fid) fprintf(fid, varargin{:}), fid,'UniformOutput',0);
end