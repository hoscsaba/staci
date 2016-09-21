function postprocess_opti(eeedummy)

    global idoszak_max
    global medence_figname nyomas_figname1 nyomas_figname2
    global AGNEVSOR PREFIX OPTI_RUN
    
    OPTI_RUN=1;  %HP 06.01.27.

    if OPTI_RUN==1
        
        logfilename='opti_log.txt';
%         data=dlmread(logfilename,'delimiter',',');
        data=dlmread(logfilename,',');  % by HP 2016.02.23.

        x=data(:,2:end-2);
        obj=data(:,end-1);

        N=ceil(sqrt(length(x(1,:))));

        subfignum=1;
        fopt=figure(10);
        for i=1:N
            for j=1:N
                if subfignum<length(x(1,:))+1
                    subplot(N,N,subfignum)
                    semilogy(x(:,subfignum),obj,'+')
                    axis([0 2.1 -inf inf])
                    xlabel(['x_',num2str(subfignum)]), grid on
                    
                    subfignum = subfignum+1;
                end
            end
        end
        print(fopt,'parameter_sensitivity','-dpdf'); orient landscape;

        [min_val,min_idx]=min(obj);
        fprintf('\n Best so far (out of %d): obj = %g\n',length(obj),min_val);
        best_x=x(min_idx,:)
        tic
        hiba_pres=foviz_objective(best_x)
        toc
    end


    % OPTI_RUN=0;       
    % medence_figname='MedenceTerfogataramok_0';
    % nyomas_figname1='Nyomasok1_0';
    % nyomas_figname2='Nyomasok2_0';
    % logfilename='opti_log_0.txt';
    % orig_x=ones(1,length(AGNEVSOR)); 
    % hiba_pres=foviz_objective(orig_x)

end