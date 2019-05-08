function  couter_of_displayed_j=plot_1D_spectra_J_data(tmp_obj,fig_num,used_color,opt)
couter_of_displayed_j=0;
figure(fig_num);hold on;
%tmp_obj=super_obj{list_1d(loop_over_spectra,1)};

if isfield(tmp_obj,'chemical_shift')
    for loop_over_peaks=1:size(tmp_obj.chemical_shift,2)
        % disp(['For Peak : ' num2str(loop_over_peaks) ' Label= ' tmp_obj.label{1,loop_over_peaks}])
        %  plot(tmp_obj.chemical_shift{1,loop_over_peaks},0,[used_color '+'])
        % text(tmp_obj.chemical_shift{1,loop_over_peaks},0,tmp_obj.label{1,loop_over_peaks},'color',used_color)
        if isfield(tmp_obj,'label')
            txt=tmp_obj.label{1,loop_over_peaks};
        else
            txt='?';
        end
        txt=[txt ' '];
        if isfield(tmp_obj,'multiplicity')
            %  if tmp_obj.multiplicity{1,loop_over_peaks}~=''
            txt=[txt num2str(tmp_obj.multiplicity{1,loop_over_peaks}) ' '];
            % end
        end
        if isfield(tmp_obj,'integral')
            if tmp_obj.integral{1,loop_over_peaks}~=0
                %  txt=[txt 'E=' num2str(tmp_obj.integral{1,loop_over_peaks}) ' '];
            end
        end
        pos_lab=loop_over_peaks;
        %%   broadened_spectrum=conv(,ones(nb_pt_for_broadening,1),'same')/nb_pt_for_broadening;
        
      
        %  disp(txt)
        list_pos=25;
        
        text(list_pos,pos_lab,txt,'color',used_color)
        %% plot j structure
        if isfield(tmp_obj,'J')
            if size(tmp_obj.J,1)>=loop_over_peaks
                curJl=zeros(1,size(tmp_obj.J,2));
                for loo_j=1: size(tmp_obj.J,2)
                    totm=tmp_obj.J{loop_over_peaks,loo_j};
                    if size(totm,2)==0
                        curJl(1,loo_j)=0;
                    else
                        curJl(1,loo_j)=totm;
                    end
                end
                curJl=sort(curJl,'descend');
                for loo_j=1: size(curJl,2)
                    curJ=curJl(1,loo_j);
                    if curJ~=0
                            plot([curJ],[ pos_lab ],[used_color 'x']);
                            text([curJ],[ pos_lab ],num2str(curJ));
                            couter_of_displayed_j=couter_of_displayed_j+1;
                    end
                end
                
            end
        end
        %% just little vertical line
            plot([0 22],[ pos_lab pos_lab ],[used_color '-']);
        
        
        
       
    end
    
    
    
end
drawnow

% print(['./plots_of_figures/M_.eps'],'-depsc');%for matlab
end
