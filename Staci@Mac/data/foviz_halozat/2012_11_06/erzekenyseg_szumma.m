function erzekenyseg_szumma

    in_fname='erzekenyseg/FV_input_2012_11_06_barlang_nelkul__1_erz.csv'

    data=read_dfile(in_fname);   
    out=cell(length(data(:,1))+2,length(data(1,:))+1);    
    
    nap_id=1;
    
    %erz=cell(1,1)
    erz.nap{1}.datum={'2012_11_06'};
    erz.nap{2}.datum={'2013_01_03'};
    erz.nap{3}.datum={'2013_07_20'};
    erz.nap{4}.datum={'2013_08_16'};
    erz.nap{5}.datum={'2013_10_01'};
    erz.nap{6}.datum={'2014_01_27'};
    erz.nap{7}.datum={'2014_02_28'};
    
%% csövek és csomópontok id-inek kigyûjtése
    for i=1:269 %sorokon léptet
       pipes(i)=data(i+1,1);
    end    
    for j=1:191
       nodes(j)=data(length(pipes)+1+j,1);
    end
    
    pipes=pipes(1:3) %For debug
    
%% napok végignézése
    for napok=1:1%7 %napok pörgetése
        for csovek=1:length(pipes) %bemenõ parmaéter: adott csõ átmérõjére vett
            %% Egy napra egy csõ értékeinek meghatározása
            
            %in_fname genrálást megírni!!!
            in_fname=['erzekenyseg/FV_input_', char(erz.nap{napok}.datum), '_barlang_nelkul__', char(pipes(csovek)), '_erz.csv'];
            data=read_dfile(in_fname)
            out=data;
            
            erz.nap{napok}.cso{csovek}.id={data(1,1)}
            
            for idoszak=1:47 %egy fájlon belül oszlopokon léptet
                q_sum=0;
                p_sum=0;
                
                for i=1:length(pipes) %sorokon léptet
                    q_sum(i)=abs(str2num(cell2mat(data(i,idoszak+1))));
                end
                sor=1+length(pipes)+length(nodes)+1;
                out(sor,idoszak+1)={sum(q_sum)};
                
                for j=1:length(nodes)
                    index=1+length(pipes)+j;
                    p_sum(j)=abs(str2num(cell2mat(data(index,idoszak+1))));
                end
                out(sor+1,idoszak+1)={sum(p_sum)};
            end
            
            out(sor,1)={'sum_q'};
            out(sor+1,1) ={'sum_p'};
            
            %% átlagok, szórások
            mean_q(napok,csovek)=mean(q_sum);
            erz.nap{napok}.cso{csovek}.q.mean=mean(q_sum);         
            erz.nap{napok}.cso{csovek}.q.std=std(q_sum);
            
            mean_p(napok,csovek)=mean(p_sum);
            erz.nap{napok}.cso{csovek}.p.mean=mean(p_sum);          
            erz.nap{napok}.cso{csovek}.p.std=std(p_sum);
            
            
            
            out_fname=[in_fname(1:end-4) '_kiertekelt.csv']
           % mentes(out_fname, out)
        end
        
        %% napi rangsorolás
%         mean_q(napok,:)
        [B rangsor_q]=sort(mean_q(napok,:))
        [B rangsor_p]=sort(mean_p(napok,:))
        
        for csovek=1:length(pipes)
           erz.nap{napok}.cso{csovek}.q.helyezes=rangsor_q(csovek)
           erz.nap{napok}.cso{csovek}.p.helyezes=rangsor_p(csovek)           
        end
        
        
        
    end
    
    
%% Egyéb
	fprintf(' kesz.\n');
end

function out = read_dfile(fname)

	fprintf('%s adatfajl olvasasa...',fname);

	fp=fopen(fname);
	line_num=1;
	while ~feof(fp)
		tmp=fgetl(fp);
		i=1;
		while ~isempty(tmp)
			[out{line_num,i},tmp]=strtok(tmp,';');
			i=i+1;
		end
		line_num=line_num+1;
	end

	fprintf(' kesz.\n');

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
