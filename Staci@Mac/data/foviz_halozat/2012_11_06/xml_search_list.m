function out_list=xml_search_list(data,type,element_list,prop)
% 	Kozvetlenul a Staci projekt fileokbol (*.spr) torteno adatolvasas

%For debug:
% 	clc
%     fname='5cso.spr';
%     type='edge';
%     element='PIPE22';
%     prop='mass_flow_rate';

%---------------------------------------------------------
%     tic
%     fname=data;
% 	data=xml_load(fname);
%   data=fname;
% fprintf('element_list hossza: %i',length(element_list));

out_list=zeros(1,length(element_list));

j=1;
switch type
    case 'edge'
        for j=1:length(element_list) %l�ptessen v�gig a keresett id-ken
            for i=1:length(data.edges(1,:)) %l�ptessen v�gig az �gakon
                
                %fprintf('data.edge: i.: %i,  %s\n',i,data.edges(1,i).edge.id);
                %fprintf('element_list: %i: %s\n\n\n *******\n',j,char(element_list{j}));
                
                if strcmp(char(element_list{j}),strtrim(data.edges(1,i).edge.id))==1  %ha az aktu�lis �g id-je megegyezik a keresettel
                    
                    switch prop
                        case 'mass_flow_rate'            %ha t�meg�ramot keres�nk
                            out = data.edges(1,i).edge.mass_flow_rate;
                        case 'volume_flow_rate'   %ha t�rfogat�ramot keres�nk
                            out = data.edges(1,i).edge.volume_flow_rate;
                        case 'diameter'
                            out = data.edges(1,i).edge.edge_spec.pipe.diameter;
                            
                        otherwise
                            error('Unknown/unimplemented property');
                    end
                    out_list(j) = str2num(out);
                    %                         i=i+1;
                    %                         if j>length(element_list)
                    %                             break
                    %                         end
                end
            end
        end
        
    case 'node'
        for j=1:length(element_list) %l�ptessen v�gig a keresett id-ken
            
            for i=1:length(data.nodes) %l�ptessen v�gig a csom�pontokon
                if strcmp(char(element_list(j)),strtrim(data.nodes(1,i).node.id)) ==1 %ha a keresett elem egyezik az aktu�lis id-vel
                    
                    switch prop
                        case 'head'       %head a kimenetre
                            out = data.nodes(1,i).node.head;
                        case 'pressure'  %nyom�s a kimenetre
                            out = data.nodes(1,i).node.pressure;
                        case 'height'    %magass�g a kimenetre
                            out = data.nodes(1,i).node.height;
                        otherwise
                            error('Unknown property');
                    end
                    out_list(j) = str2double(out);
                    %                 j=j+1;
                    %                 if j>length(element_list)
                    %                     break
                    %                 end
                end
            end
        end
        
        
    otherwise
        error('Unknown type');
end
%     toc
%     out_list(j)=str2num(out);

end