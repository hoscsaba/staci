function hiba_pres=foviz_objective(x)

    global logfname spr DEBUG AGNEVSOR OPTI_RUN
    global staci_cmd_string
    global opti_logfile_name fun_eval
    global idoszak_max
    global medence_figname nyomas_figname1 nyomas_figname2
    global FOVIZ_DFILE_NAME fname

    staci_cmd_string = '/home/hoscsaba/Dropbox/staci/kod/new_staci/Staci@Wylie/new_staci';
%staci_cmd_string = '/Users/hoscsaba/Dropbox/staci/kod/new_staci/Staci@Mac/new_staci';
%---------------------

CHANGE_DFILES   = 1;
RUN_DFILES     = 1;
POSTPROCESSING_POOL = 1;
POSTPROCESSING_PRESSURE = 1;
DRAW_FIGS=1;

DEBUG=1;

%---------------------

spr='alapzonamodell9.spr';
fname = FOVIZ_DFILE_NAME;
logfname = 'log.txt';
fp=fopen(logfname,'w');
fclose(fp);

data = read_dfile(fname);
data = sort_data(data);


if CHANGE_DFILES==1

    fp=fopen(logfname,'a');
    mfprintf([1 fp],'\n\n   CHANGING PIPE DIAMETERS...\n======================================');
    fclose(fp);
    
% %     for i=1:length(AGNEVSOR)
% % 
% %         spr=[fname(1:end-4),'_1.spr'];
% %         D=get_pipe_diameter(spr,AGNEVSOR{i});
% %         Dnew = D*x(i);
% %         
% %         fp=fopen(logfname,'a');
% %         mfprintf([1 fp],'\n\t Pipe %s: %g -> %g',AGNEVSOR{i},D,Dnew);
% %         
% %         for idoszak=1:idoszak_max
% %             spr=[fname(1:end-4),'_',num2str(idoszak),'.spr'];
% %             spr_new=[fname(1:end-4),'_',num2str(idoszak),'_mod.spr'];
% %             if exist(spr_new,'file')
% %                 delete(spr_new);
% %             end
% %             set_pipe_diameter(spr,AGNEVSOR{i},Dnew,spr_new);
% %         end
% %         fclose(fp);
% %     end

% Majdnem jó volt a fájlokban átmérőket változtató eljárás, de csak az
% utolsóként beálltott csővet változtatta meg.
    
    for idoszak=1:idoszak_max  %Régi fájlok törlése, újak létrehozása
        spr=[fname(1:end-4),'_',num2str(idoszak),'.spr'];
        spr_new=[fname(1:end-4),'_',num2str(idoszak),'_mod.spr'];
            if exist(spr_new,'file')
                delete(spr_new);
            end
            D=get_pipe_diameter(spr,AGNEVSOR{1});
            set_pipe_diameter(spr,AGNEVSOR{1},D,spr_new); %csak fájl létrehozáshoz!!
    end

    for i=1:length(AGNEVSOR)%átmérő átállítás

        spr=[fname(1:end-4),'_1.spr'];
        D=get_pipe_diameter(spr,AGNEVSOR{i});
        Dnew = D*x(i);
%         pause
        fp=fopen(logfname,'a');
        mfprintf([1 fp],'\n\t Pipe %s: %g -> %g',AGNEVSOR{i},D,Dnew);
        
        fclose(fp);
        for idoszak=1:idoszak_max
%             spr=[fname(1:end-4),'_',num2str(idoszak),'.spr'];
            spr_new=[fname(1:end-4),'_',num2str(idoszak),'_mod.spr'];
            set_pipe_diameter(spr_new,AGNEVSOR{i},Dnew,'temp.spr');
            
            delete(spr_new);
            movefile('temp.spr',spr_new);
        end
        
    end
    

%     for i=1:length(AGNEVSOR) %check 
%         for idoszak=1:idoszak_max
%                    
%             spr=[fname(1:end-4),'_1.spr'];
%             D=get_pipe_diameter(spr,AGNEVSOR{i});
%             Dnew = D*x(i);
%         
%             spr_check=[fname(1:end-4),'_',num2str(idoszak),'_mod.spr'];
%             Dcheck = get_pipe_diameter(spr_check,AGNEVSOR{i});
%             mfprintf([1 fp],'\n %g, %s: Dnew=%g, Dcheck=%g',idoszak,AGNEVSOR{i},Dnew,Dcheck);
% %             pause
%             
%         end
% %         pause
%     end
%     fclose(fp);
end
% pause

if RUN_DFILES==1

    fp=fopen(logfname,'a');
    mfprintf([1 fp],'\n\n   RUNNING SIMULATIONS...\n=================================');
    fclose(fp);
    
    data = read_dfile(fname);
    data = sort_data(data);
    
    for idoszak=1:idoszak_max
        fp=fopen(logfname,'a');
        mfprintf([1 fp],'\n\n   IDOSZAK: %g\n-------------------------------',idoszak);
        fclose(fp);
        
        % Filenev beallitasa
        spr_new=[fname(1:end-4),'_',num2str(idoszak),'_mod.spr'];
        
        % Staci futtatasa
        is_ok=run_staci(spr_new);
        if (is_ok~=1)==1
            error('SIKERTELEN SZAMITAS!!!');
        end
        
        % Szamitott medenceterfogataramok gyujtese
        %get_pool_flow_rates(idoszak,data,spr_new);
    end
end

if POSTPROCESSING_POOL==1

    data = read_dfile(fname);
    data = sort_data(data);
    
    fp=fopen(logfname,'a');
    mfprintf([1 fp],'\n\n   POSTPROCESSING...\n=================================');
    fclose(fp);
    
    
    % MEDENCE TERFOGATARAMOK
    
    for idoszak=1:idoszak_max
        fp=fopen(logfname,'a');
        mfprintf([1 fp],'\n\n   IDOSZAK: %g\n-------------------------------',idoszak);
        fclose(fp);
        
        % Filenev beallitasa
        spr_new=[fname(1:end-4),'_',num2str(idoszak),'_mod.spr'];
        
        % Szamitott medenceterfogataramok gyujtese
        [QpoolHDS(idoszak,:),QpoolFVM(idoszak,:)]=get_pool_flow_rates(idoszak,data,spr_new);
    end
    
    tveclength=length(QpoolFVM(:,1));
    tvec=(1:1:tveclength)/2;
    f1=figure(1);
    for i=1:length(QpoolHDS(1,:))

        subplot(4,2,i)
        minQ=min(min(QpoolHDS(:,i)),min(QpoolFVM(:,i)))/1000;
        maxQ=max(max(QpoolHDS(:,i)),max(QpoolFVM(:,i)))/1000;
        plot(tvec,QpoolHDS(:,i)/1000,'r+-',tvec,QpoolFVM(:,i)/1000,'kx--');
        if i==1
            title([data.med{i}.name,'+: HDS, x: FV'],'Interpreter','none')
        else
            title(data.med{i}.name,'Interpreter','none')
        end
        grid on, xlabel('t [h]'),ylabel('Q [ezer m3/h]');
        %legend('HDS','FV')
        ylim([abs(minQ)*1.2*sign(minQ),abs(maxQ)*1.2*sign(maxQ)])
        xlim([0 24])
    end
    
    %print(f1,medence_figname,'-dpng','-r800');
%     print(f1,medence_figname,'-dpdf'); orient landscape;
    medence_figname_2=['MedenceTerfogataramok_',num2str(fun_eval),'.pdf'];
    print(f1,medence_figname_2,'-dpdf'); orient landscape;
end


% NYOMASOK
hiba_pres=0;
if POSTPROCESSING_PRESSURE==1

    % NYOMASOK
    hiba_pres=0;
    for idoszak=1:idoszak_max

        % Filenev beallitasa
        spr_new=[fname(1:end-4),'_',num2str(idoszak),'_mod.spr'];
        
        % Szamitott medenceterfogataramok gyujtese
        [PresHDS(idoszak,:),PresFVM(idoszak,:),nodenames]=get_spec_pres(idoszak,data,spr_new);
        
%         x1=PresHDS(idoszak,:);
%         x2=PresFVM(idoszak,:);
        %x1=x1-mean(x1);
        %x2=x2-mean(x2);
%         
%         hiba_pres=hiba_pres+norm(abs(x1-x2));
%         fp=fopen(logfname,'a');
%         mfprintf([1 fp],'\n\n   hiba_pres = %g\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~',hiba_pres);
%         fclose(fp);
    end
    
    %%by HP hiba_pres átírása
    for pont=1:length(PresHDS(1,:))
        %meanHDS=mean(PresHDS(:,pont));
        %meanFVM=mean(PresFVM(:,pont));
        
        PresHDS_tr=PresHDS(1:(end-1),pont);
        PresFVM_tr=PresFVM(1:(end-1),pont);
%         length(PresHDS)
        
%         pause
%         pause
        
        meanHDS=mean(PresHDS_tr);
        meanFVM=mean(PresFVM_tr);
        
        nullHDS=PresHDS_tr-meanHDS;
        nullFVM=PresFVM_tr-meanFVM;
        
        kul=nullHDS-nullFVM;
        
        hiba=norm(kul);
        hiba_pres=hiba_pres+hiba;
%         pause
        
        fp=fopen(logfname,'a');
        mfprintf([1 fp],'\n\n   hiba_pres = %g\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~',hiba_pres);
        fclose(fp);
    end
    
    %% end of hiba_pres átírása
    
    if DRAW_FIGS==1
        %tveclength=length(QpoolFVM(:,1));
        %tvec=(1:1:tveclength)/2;
        f2=figure(2);orient landscape;
        f3=figure(3);orient landscape;
        for i=1:length(PresHDS(1,:))
            if i<11
                figure(2)
                subplot(5,2,i)
            else
                figure(3)
                subplot(5,2,i-10)
            end
            minPres=min(min(PresHDS(:,i)),min(PresFVM(:,i)));
            maxPres=max(max(PresHDS(:,i)),max(PresFVM(:,i)));
            plot(tvec,PresHDS(:,i),'r+-',tvec,PresFVM(:,i),'kx--');
            title(nodenames{i},'Interpreter','none'), grid on, xlabel('t [h]'),ylabel('H [m]');
            %legend('HDS','FV')
            if minPres<0
                mul=0.9;
            else
                mul=1.2;
            end
            %ylim([abs(minPres)*mul*sign(minPres),abs(maxPres)*1.2*sign(maxPres)])
            xlim([0 24])
        end
        %print(f2,nyomas_figname1,'-dpng','-r800');
        %print(f3,nyomas_figname2,'-dpng','-r800');
        
%         print(f2,nyomas_figname1,'-dpdf');
%         print(f3,nyomas_figname2,'-dpdf');

        nyomas_figname1_2=['Nyomasok1_',num2str(fun_eval),'.pdf'];
        nyomas_figname2_2=['Nyomasok2_',num2str(fun_eval),'.pdf'];
        print(f2,nyomas_figname1_2,'-dpdf');
        print(f3,nyomas_figname2_2,'-dpdf');
    end
end

if OPTI_RUN==1

    fun_eval=fun_eval+1;

    fp=fopen(opti_logfile_name,'a');
    fprintf(fp,'%2g,',fun_eval);
    for i=1:length(x)
        fprintf(fp,'%6.4e, ',x(i));
    end
    fprintf(fp,'%6.4e, %s\n',hiba_pres,num2str(clock));
    fclose(fp);
end

delete('*.spr.ros','*.spr.rps','*.spr.rrs');

end

%------------------------------------------------
function set_pipe_diameter(fname,edgename,value,fname_new)
    global staci_cmd_string
    staci_cmd = sprintf('system(''%s -m %s -e %s -p diameter -n %g -o %s'')',...
        staci_cmd_string,fname,edgename,value,fname_new);
    res=strtrim(evalc(staci_cmd));

%copyfile('tmp.spr',fname);
%delete('tmp.spr');
end


%------------------------------------------------
function out=get_pipe_diameter(fname,edgename)
    global staci_cmd_string
    staci_cmd = sprintf('system(''%s -g %s -e %s -p diameter'')',staci_cmd_string,fname,edgename);

    res=strtrim(evalc(staci_cmd));

    idx1=strfind(res,'value');
    idx2=strfind(res,'Timing');
    [tok,rem]=strtok(res(idx1+6:idx2-1),':');
    out=str2num(strtrim(rem(2:end)));
end


%------------------------------------------------
function [H_HDS,H_FV,nodenames]=get_spec_pres(idoszak,data,spr_new)

    global logfname spr
    global FOVIZ_DFILE_NAME

    fp=fopen(logfname,'a');
    mfprintf([1 fp],'\n\n   Nyomasok ellenorzese: (idoszak=%g) ',idoszak);

    data=read_dfile(['Eredmenyek_',FOVIZ_DFILE_NAME]);
    save('tmp.mat','data');
    kk=1;
        
    for i=2:length(data)
        FV_ID=data{i}{1};
        Staci_ID=data{i}{2};
             FV_ID_v(kk)=cell(1); %by hp
             
        if length(Staci_ID)>2
            
            FV_ID_v(kk)={FV_ID}; %by hp            
%             get_node_head(spr_new,Staci_ID);  %% commented by HP 03.07.
            H_FV(kk) = str2num(data{i}{4+idoszak});
%             H_HDS(kk) = get_node_head(spr_new,Staci_ID) %commented by HP

            nodenames{kk}=Staci_ID;
%             mfprintf([1 fp],'\n\t FV: %s, Staci alias: %s, nyomas: FV : %g vom, HDS : %g vom',FV_ID,Staci_ID,H_FV(kk),H_HDS(kk)); 
            kk=kk+1;
        end
    end
        
    data_xml=xml_load(spr_new); %by hp
    H_HDS=xml_search_list(data_xml,'node',nodenames,'head'); %by hp
    
    for k=1:length(H_HDS)
        mfprintf([1 fp],'\n\t FV: %s, Staci alias: %s, nyomas: FV : %g vom, HDS : %g vom',char(FV_ID_v{k}),char(nodenames{k}),H_FV(k),H_HDS(k));
    end
%     pause %by hp 03.07.
    
    fclose(fp);

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
%  global logfname spr
% 
%     fp=fopen(logfname,'a');
%     mfprintf([1 fp],'\n\n   Medencek staci terfogatarama:');
% 
%     for i=1:length(data.med)
%         mfprintf([1 fp],'\n     %s :',add_whitespace(data.med{i}.name,14));
%         QpoolHDS(i)=get_edge_flowrate(fname,data.med{i}.staci_name);
%         QpoolFVM(i)=data.med{i}.A*(data.med{i}.level(idoszak+1)-data.med{i}.level(idoszak))*2;
%         mfprintf([1 fp],'  %+6.2f m3/h, Foviz adat: %+6.2f',QpoolHDS(i),QpoolFVM(i));
%         mfprintf([1 fp],'\t  %+6.2f = ( %g - %g ) * %g * 2',QpoolFVM(i),...
%             data.med{i}.level(idoszak+1),data.med{i}.level(idoszak),data.med{i}.A);
% 
%     end
%     fclose(fp);
%     pause

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
    global staci_cmd_string

    fp=fopen(logfname,'a');
    fname_rrs=[fname,'.rrs'];
    if exist(fname_rrs,'file')==2
        mfprintf([1 fp],'\n\n Deleting %s ...',fname_rrs);
        delete(fname_rrs);
    end

    staci_cmd = sprintf('system(''%s -s %s'')',staci_cmd_string,fname);
    mfprintf([1 fp],'\n\n Running %s ...',staci_cmd);

    fclose(fp);

    res=strtrim(evalc(staci_cmd));

    fp=fopen(logfname,'a');
    mfprintf([1 fp],'%s',res);
    fclose(fp);

    is_ok=strcmp(strtrim(fileread(fname_rrs)),'ok');
end

%------------------------------------------------
function set_pool_water_level(fname,poolname,value)
    global staci_cmd_string

    staci_cmd = sprintf('system(''%s -m %s -e %s -p water_level -n %g -o %s'')',...
        staci_cmd_string,fname,poolname,value,'tmp.spr');
    res=strtrim(evalc(staci_cmd));

    copyfile('tmp.spr',fname);
    delete('tmp.spr');
end

%------------------------------------------------
function set_valve_opening(fname,valvename,value)
    global staci_cmd_string
    staci_cmd = sprintf('system(''%s -m %s -e %s -p position -n %g -o %s'')',...
        staci_cmd_string,fname,valvename,value,'tmp.spr');
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

set_zone_consumptions(fname,data,'A_',Q_A_kivant/Q_A);
set_zone_consumptions(fname,data,'K_',Q_K_kivant/Q_K);

% ELLENORZES
%Q_A_uj=get_zone_consumptions(fname,data,'A_');
%Q_K_uj=get_zone_consumptions(fname,data,'K_');

% SPECIALIS CSOMOPONTOK
fp=fopen(logfname,'a');
mfprintf([1 fp],'\n   Specialis csompontok:');
fclose(fp);

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
            fp=fopen(logfname,'a');
            mfprintf([1 fp],'\n \t %s:\t %6.2f m3/h (betap)',node_name,Qcsp);
            fclose(fp);
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
        Qout=get_node_consumption(fname,data.tovabbemeles{i}.name)*data.tovabbemeles{i}.mul;
        Qzona=Qzona-Qout;
        node_name = add_whitespace(data.tovabbemeles{i}.name,22);
        mfprintf([1 fp],'\n \t %s:\t %6.1f m3/h (tovabbemeles), Qzona=%6.1f m3/h',node_name,Qout,Qzona);
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
            if DEBUG>1
                mfprintf([1 fp],'\n \t %s:\t %6.2f m3/h',node_name,Qcsp);
            end
            if sum(strcmp(node_list{i},data.spec_nodes))==0
                Q=Qcsp*mul;
                set_node_consumption(fname,node_list{i},Q);
                if DEBUG>1
                    mfprintf([1 fp],'\t -> %6.2f m3/h',Q);
                    Q1=get_node_consumption(fname,node_name);
                    mfprintf([1 fp],' (ellenorzes: Q=%g m3/h)',Q1);
                end
            else
                if DEBUG>1
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
    global staci_cmd_string
    staci_cmd = sprintf('system(''%s -m %s -e %s -p demand -n %g -o %s'')',...
        staci_cmd_string,fname,nodename,value,'tmp.spr');
    res=strtrim(evalc(staci_cmd));

    copyfile('tmp.spr',fname);
    delete('tmp.spr');
end

%------------------------------------------------
function Qcsp=get_node_consumption(fname,nodename)
    global staci_cmd_string
    staci_cmd = sprintf('system(''%s -g %s -e %s -p demand'')',staci_cmd_string,fname,nodename);

    res=strtrim(evalc(staci_cmd));

    idx1=strfind(res,'value');
    idx2=strfind(res,'Timing');
    [tok,rem]=strtok(res(idx1+6:idx2-1),':');
    Qcsp=str2num(strtrim(rem(2:end)));
end

%------------------------------------------------
function pcsp=get_node_head(fname,nodename)
    global staci_cmd_string

    staci_cmd = sprintf('system(''%s -g %s -e %s -p head'')',staci_cmd_string,fname,nodename);

    res=strtrim(evalc(staci_cmd));

    idx1=strfind(res,'value');
    idx2=strfind(res,'Timing');
    [tok,rem]=strtok(res(idx1+6:idx2-1),':');
    pcsp=str2num(strtrim(rem(2:end)));
end


%------------------------------------------------
function Qag=get_edge_flowrate(fname,edgename)
    global staci_cmd_string
    staci_cmd = sprintf('system(''%s -g %s -e %s -p mass_flow_rate'')',staci_cmd_string,fname,edgename);

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
            if DEBUG>1
                mfprintf([1 fp],'\n \t %s:\t %6.2f m3/h',node_name,Qcsp);
            end

            if sum(strcmp(node_list{i},data.spec_nodes))==0
                Q=Q+Qcsp;
                if DEBUG>1
                    mfprintf([1 fp],'\t zona : %6.2f m3/h',Q);
                end
            else
                if DEBUG>1
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

function mfprintf(fid, varargin)
    arrayfun(@(fid) fprintf(fid, varargin{:}), fid,'UniformOutput',0);
end