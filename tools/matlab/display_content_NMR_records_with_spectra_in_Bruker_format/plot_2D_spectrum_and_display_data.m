function  plot_2D_spectrum_and_display_data(super_obj,object_num,used_color,opt)
figure(object_num);hold on;

tmp_obj=super_obj{object_num};


%%%%%%%%%%

min_chem_shift_spectrum1=10000000;
max_chem_shift_spectrum1=-10000000000;
min_chem_shift_spectrum2=10000000;
max_chem_shift_spectrum2=-10000000000;
if isfield(tmp_obj,'signalf1')
    for loop_over_peaks=1:size(tmp_obj.signalf1,2)
        chem_shift1=0;
        chem_shift2=0;
        % disp(['For Peak : ' num2str(loop_over_peaks) ' Label= ' tmp_obj.label{1,loop_over_peaks}])
        %  plot(tmp_obj.chemical_shift{1,loop_over_peaks},0,[used_color '+'])
        % text(tmp_obj.chemical_shift{1,loop_over_peaks},0,tmp_obj.label{1,loop_over_peaks},'color',used_color)
        txt='';
        % analyse relevance
        lab1=tmp_obj.signalf1{1,loop_over_peaks};
        lab2=tmp_obj.signalf2{1,loop_over_peaks};
        if (size(lab1,2)>0) && (size(lab2,2)>0)
            index1=get_index_of_label(super_obj{1,1},lab1);
            index2=get_index_of_label(super_obj{1,1},lab2);
            chem_shift1=super_obj{1,1}.chemical_shift{1,index1};
            chem_shift2=super_obj{1,1}.chemical_shift{1,index2};
            
            aaa1= super_obj{1}.structure.ref_from_label_number_to_nb_bond_table(index1);
            aaa2= super_obj{1}.structure.ref_from_label_number_to_nb_bond_table(index2);
            li_nb=super_obj{1}.structure.nb_bond_between_atoms_including_implicit_H(:,aaa1,aaa2);
            %
            %  txt=[txt tmp_obj.signalf1{1,loop_over_peaks} '/'];
            % txt=[txt tmp_obj.signalf2{1,loop_over_peaks}];
            txt=[txt lab1 '/'];
            txt=[txt lab2];
            for loi=1:size(li_nb,1)
                if li_nb(loi)
                    txt=[txt ' ' num2str(loi) 'J'];
                end
            end
            
            %                                 for loop_over_assigned_signals=1:size(super_obj{1,1}.label_signal,2)
            %                                     if size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},2)>0
            %                                         if size(tmp_obj.signalf1{1,loop_over_peaks},2)>0
            %                                             %  disp(['compare '      super_obj{1,1}.label_signal{1,loop_over_assigned_signals} ' '    tmp_obj.signalf1{1,loop_over_peaks}])
            %                                             %  disp(['compare '      num2str(size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},1)) ' '     num2str(size(tmp_obj.signalf1{1,loop_over_peaks},1))])
            %                                             %  disp(['compare '       num2str(size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},2)) ' '     num2str(size(tmp_obj.signalf1{1,loop_over_peaks},2))])
            %                                             if strcmp(super_obj{1,1}.label_signal{1,loop_over_assigned_signals} , tmp_obj.signalf1{1,loop_over_peaks})
            %                                                 list_pos1=super_obj{1,1}.chemical_shift{1,loop_over_assigned_signals};
            %                                             end
            %                                         end
            %                                     end
            %                                 end
            %                                 for loop_over_assigned_signals=1:size(super_obj{1,1}.label_signal,2)
            %                                     if size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},2)>0
            %                                         if size(tmp_obj.signalf2{1,loop_over_peaks},2)>0
            %                                             %  disp(['compare '      super_obj{1,1}.label_signal{1,loop_over_assigned_signals} ' '    tmp_obj.signalf1{1,loop_over_peaks}])
            %                                             %  disp(['compare '      num2str(size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},1)) ' '     num2str(size(tmp_obj.signalf1{1,loop_over_peaks},1))])
            %                                             %  disp(['compare '       num2str(size(super_obj{1,1}.label_signal{1,loop_over_assigned_signals},2)) ' '     num2str(size(tmp_obj.signalf1{1,loop_over_peaks},2))])
            %                                             if strcmp(super_obj{1,1}.label_signal{1,loop_over_assigned_signals} , tmp_obj.signalf2{1,loop_over_peaks})
            %                                                 list_pos2=super_obj{1,1}.chemical_shift{1,loop_over_assigned_signals};
            %                                             end
            %                                         end
            %                                     end
            %                                 end
        end
        if (chem_shift2~=0) && (chem_shift1~=0)
            plot(chem_shift2,chem_shift1,[used_color '+'])
            text(chem_shift2,chem_shift1,txt,'color',used_color)
            
            
            if opt.validate(1,4)>0% opt.validate 2D correlations
                
                xlim(chem_shift2+[-0.25 0.25]);
                ylim(chem_shift1+[-0.25 0.25]);
                
                choice = questdlg(['Is the signal assigned to ' txt ' really present?'], 'Validation of 2D correlations', 'Yes','No','No, add comment','Yes');
                opt.validate(1,4)=opt.validate(1,4)-1;
            end
            
            
            min_chem_shift_spectrum1=min(min([min_chem_shift_spectrum1 chem_shift1]));
            max_chem_shift_spectrum1=max(max([max_chem_shift_spectrum1 chem_shift1]));
            min_chem_shift_spectrum2=min(min([min_chem_shift_spectrum2 chem_shift2]));
            max_chem_shift_spectrum2=max(max([max_chem_shift_spectrum2 chem_shift2]));
        end
        
        
        
    end
    %set boudaries of axis
    margin1=0.2;fp1=power(10,-1);
    margin2=0.2;fp2=power(10,-1);
    if abs(min_chem_shift_spectrum1-max_chem_shift_spectrum1)>10
        margin1=1;fp1=power(10,0);
    end
    if abs(min_chem_shift_spectrum1-max_chem_shift_spectrum1)>100
        margin1=10;fp1=power(10,1);
    end
    if abs(min_chem_shift_spectrum2-max_chem_shift_spectrum2)>10
        margin2=1;fp2=power(10,0);
    end
    if abs(min_chem_shift_spectrum2-max_chem_shift_spectrum2)>100
        margin2=10;fp2=power(10,1);
    end
    min_chem_shift_spectrum1=min_chem_shift_spectrum1 -margin1;min_chem_shift_spectrum1=round(min_chem_shift_spectrum1/fp1)*fp1;
    max_chem_shift_spectrum1=max_chem_shift_spectrum1 +margin1;max_chem_shift_spectrum1=round(max_chem_shift_spectrum1/fp1)*fp1;
    min_chem_shift_spectrum2=min_chem_shift_spectrum2 -margin2;min_chem_shift_spectrum2=round(min_chem_shift_spectrum2/fp2)*fp2;
    max_chem_shift_spectrum2=max_chem_shift_spectrum2 +margin2;max_chem_shift_spectrum2=round(max_chem_shift_spectrum2/fp2)*fp2;
    xlim([min_chem_shift_spectrum2 max_chem_shift_spectrum2])
    ylim([min_chem_shift_spectrum1 max_chem_shift_spectrum1])
end
drawnow

end