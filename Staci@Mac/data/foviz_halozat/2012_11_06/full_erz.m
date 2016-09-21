function full_erz

%% Init
	path1 = getenv('PATH');
	path1 = [path1 '/home/hoscsaba/Dropbox/staci/kod/new_staci/Staci@Wylie'];
	setenv('PATH', path1);

    spr_file='FV_input_2012_11_06_barlang_nelkul_1.spr';
%     spr_file='5csoxyxyxyxyxyxy.spr';
%     futtat(spr_file);
    pipes=staci_list(spr_file, 'Cso');
    elements=staci_list_all(spr_file);
    
    ered=cell(length(elements)+1,48);
    
    for k=1:length(elements) %elsõ oszlop feltöltése
        ered(1+k,1)=elements(k);
    end
   
    for i=1:length(pipes) %egy fájlba menteni
        fprintf('%s\n',char(pipes(i)))
        ered(1,1)=pipes(i);
        
        %% érintõk gyûjtése
        
        for idoszak=1:47 %egy oszlop 47
            ered(1,idoszak+1)={idoszak};
            spr_akt=[spr_file(1:end-5),num2str(idoszak),'.spr'];
%           
            staci_erz(spr_akt,char(pipes(i)));
            out=import_dxdmu(elements);
            
            for k=1:length(elements) %egy sor
                ered(k+1,idoszak+1)=out(k,2);
            end
        end
        
        ered_file=[spr_file(1:end-5),'_' char(pipes{i}) ,'_erz.csv'];         
        mentes(ered_file,ered);
        movefile(ered_file,'../2012_11_06/erzekenyseg');
    end


end

function mentes(fname,adatok)
fid = fopen(fname,'w');
delimiter=';';
    for z=1:size(adatok,1)
        for s=1:size(adatok,2)
            
            var = eval(['adatok{z,s}']);
            if size(var,1) == 0
                var = '';
            end
            if isnumeric(var) == 1
                var = num2str(var);
            end
            
            fprintf(fid,var);            
            if s ~= size(adatok,2)
                fprintf(fid,[delimiter]);
            end
        end
        fprintf(fid,'\n');
    end
fclose(fid);
end

function staci_erz(fname,element)    
global staci_cmd_string
    staci_cmd=sprintf('system(''new_staci -r %s -e %s -p diameter'')',fname,element);
%     fprintf('staci parancs: %s', staci_cmd);
    res=strtrim(evalc(staci_cmd));
end

function out=import_dxdmu(elements)

    out=cell(length(elements),2);
	fp=fopen('dxdmu.txt');
    kuka=fgetl(fp); %fejléc kukázása
    
    while ~feof(fp)
       line=fgetl(fp);
       index_e=strfind(line,'@')+2;
       index_v=strfind(line,'); ')-1;
       elem=line(index_e:index_v);
%        line(index_v+4:index_v+13)
       erinto=str2double(line(index_v+4:index_v+13));
       for j=1:length(elements)
           if strcmp(elem,char(elements{j}))==1                        
                out(j,2)={erinto};
           end
       end
    end    

%     for i=1:length(elements)
%        out(i,1)={elem};
%        fprintf('%s %s \n',char(out{i,1}),out{i,2}); 
%     end
    fclose(fp);
end

function out=import_dxdmu_R(elements)

	fp=fopen('dxdmu.txt');
	line_num=1;
    index=[];
      
    for pp=1:length(elements) %sorokon léptet
       
         keres=[char(elements(pp)),');'];
         out(line_num,1)=elements(pp);
          
        while isempty(index)
            tmp=fgetl(fp);
            index=strfind(tmp,keres)+length(keres);            
        end        
        
        tmp_str=tmp(index:end);
        index_v=strfind(tmp_str,';')+index;
        out(line_num,2)=strtrim({tmp( index:index_v )});
%         pause
        line_num=line_num+1;
        index=[];
    end

    for i=1:(length(nodes)+length(pipes))
       fprintf('%s %s \n',char(out{i,1}),out{i,2}); 
    end
end

function [element_list] = staci_list_all(fname) %2. argumentum element volt

    staci_l=sprintf('system(''new_staci -l %s'')',fname);
    list=evalc(staci_l)
          
     %visszaérkezõ adatok trimmelése
     ind_eleje=strfind(list,'deleted.')+length('deleted.')
     ind_vege=strfind(list,'Timing')-1
%      fprintf('vege: %s, vege+1: %s',list(ind_vege), list(ind_vege+1));
     list=list(ind_eleje:ind_vege)
     
     rem=list(2:end);
     i=1;
     while length(rem)>1
        [type_act, rem]=strtok(rem,';');
        rem=rem(2:end);
        [elem, rem]=strtok(rem,';');
%         if strcmp(type,type_act)==1
           element_list{i}=strtrim(elem)
           i=i+1;
           %fprintf('elem: %s\n',elem);
%         end
        rem=rem(3:end);
     end
end

function futtat(fname)    
    staci_cmd=sprintf('system(''new_staci -s %s'')',fname);
    res=strtrim(evalc(staci_cmd));
%         pause
end

function [element_list] = staci_list(fname, type) %2. argumentum element volt

    staci_l=sprintf('system(''new_staci -l %s'')',fname);
    list=evalc(staci_l);
          
     %visszaérkezõ adatok trimmelése
     ind_eleje=strfind(list,'deleted.')+length('deleted.')
     ind_vege=strfind(list,'Timing')-1
%      fprintf('vege: %s, vege+1: %s',list(ind_vege), list(ind_vege+1));
     list=list(ind_eleje:ind_vege);
     
     rem=list(2:end);
     i=1;
     while length(rem)>1
        [type_act, rem]=strtok(rem,';');
        rem=rem(2:end);
        [elem, rem]=strtok(rem,';');
        if strcmp(type,type_act)==1
           element_list{i}=strtrim(elem);
           i=i+1;
           %fprintf('elem: %s\n',elem);
        end
        rem=rem(3:end);
     end
end
%------------------------------------------------
function set_element_parameter(fname,elementname,parameter,value)
	staci_cmd = sprintf('system(''new_staci -m %s -e %s -p %s -n %g -o %s'')',...
		fname,elementname,parameter,value,'tmp.spr');
	res=strtrim(evalc(staci_cmd));

	copyfile('tmp.spr',fname);
	delete('tmp.spr');
end

%------------------------------------------------
function value=get_element_parameter(fname,elementname,parameter)
	staci_cmd = sprintf('system(''new_staci -g %s -e %s -p %s'')'...
        ,fname,elementname,parameter);

	res=strtrim(evalc(staci_cmd));

	idx1=strfind(res,'value');
	idx2=strfind(res,'Timing');
	[tok,rem]=strtok(res(idx1+6:idx2-1),':');
	value=str2num(strtrim(rem(2:end)));
end

%------------------------------------------------
function out = read_dfile(fname)

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


function printp(h,figname)
      set(h,'PaperSize',fliplr(get(h,'PaperSize')));
      set(h,'PaperUnits','centimeters');
      set(h,'Units','centimeters');
      pos=get(h,'Position');
      set(h,'PaperSize',[pos(3) pos(4)]);
      set(h,'PaperPosition',[0 0 pos(3) pos(4)]);
      print(h,'-dpdf','-r200',figname);    
%       print(h,'-djpeg','-r200',figname);
end