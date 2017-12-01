function  plot_1D_spectrum_and_display_data(tmp_obj,fig_num,used_color,opt)

figure(fig_num);hold on;
%tmp_obj=super_obj{list_1d(loop_over_spectra,1)};
min_chem_shift_spectrum=10000000;
max_chem_shift_spectrum=-10000000000;
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
        pos_lab=max(max(tmp_obj.spectrum.spectrum))/2;
        nb_pt_for_broadening=1+round(0.025/(max(max(tmp_obj.spectrum.scale2))-min(min(tmp_obj.spectrum.scale2)))*tmp_obj.spectrum.si2);
        %%   broadened_spectrum=conv(,ones(nb_pt_for_broadening,1),'same')/nb_pt_for_broadening;
        [del, index]=min(abs(tmp_obj.spectrum.scale2-   tmp_obj.chemical_shift{1,loop_over_peaks}));
        from=index-nb_pt_for_broadening;
        to=index+nb_pt_for_broadening;
        if from<1, from=1;end
        if to>tmp_obj.spectrum.si2, to=si2;end
        pos_lab=max(max(tmp_obj.spectrum.spectrum(from:to,1)));
        
        if isfield(tmp_obj,'intensity')
            if size(tmp_obj.intensity,2)>=loop_over_peaks
                if tmp_obj.intensity{1,loop_over_peaks}~=0
                    pos_lab=tmp_obj.intensity{1,loop_over_peaks};
                end
            end
        end
        %  disp(txt)
        list_pos=tmp_obj.chemical_shift{1,loop_over_peaks};
        
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
                        shiftJ=0.5*curJ/tmp_obj.spectrum.sf2;
                        alist_pos=[list_pos list_pos];
                        list_pos=[list_pos+shiftJ list_pos-shiftJ];
                        n_pos_lab=pos_lab-curJ*pos_lab*0.005;
                        for lol=1:size(alist_pos,2)
                            plot([alist_pos(1,lol) list_pos(1,lol)],[ pos_lab n_pos_lab],[used_color '-']);
                        end
                        pos_lab=n_pos_lab;
                    end
                end
                
            end
        end
        %% just little vertical line
        alist_pos=[list_pos];
        n_pos_lab=pos_lab-4*pos_lab*0.005;
        for lol=1:size(alist_pos,2)
            plot([alist_pos(1,lol) list_pos(1,lol)],[ pos_lab n_pos_lab],[used_color '-']);
        end
        
        
        
        min_chem_shift_spectrum=min(min([min_chem_shift_spectrum list_pos]));
        max_chem_shift_spectrum=max(max([max_chem_shift_spectrum list_pos]));
    end
    
    
    
    %% ask to opt.validate each multiplet...
    for loop_over_peaks=1:size(tmp_obj.chemical_shift,2)
        if isfield(tmp_obj,'multiplicity')
            if opt.validate(1,3)>0
                
                list_pos=tmp_obj.chemical_shift{1,loop_over_peaks};
                %  if tmp_obj.multiplicity{1,loop_over_peaks}~=''
                txt=[ num2str(tmp_obj.multiplicity{1,loop_over_peaks}) ' '];
                % end
                
                figure(list_1d(loop_over_spectra,1));hold on;
                xlim(list_pos+[-0.05 0.05]);
                
                %%%    plot(list_pos,0*list_pos,'k.')%just to bring figure to top in octove
                %%%  drawnow('expose')
                %%%   refresh
                %%%    title('test if this brings fig to top...')
                labb='';
                if isfield(tmp_obj,'label')
                    labb=tmp_obj.label{1,loop_over_peaks};
                end
                % see if couplings are available for this...
                coupling_present=0;
                if size(tmp_obj.J,1)>=loop_over_peaks
                    if  size(tmp_obj.J{loop_over_peaks,1}) >0
                        coupling_present=1;
                    end
                end
                if coupling_present
                    choice = questdlg(['Is this multiplet of ' labb ' really a ' txt ' with these splittings ?'], 'Validation of 1D multiplet', 'Yes','No','No, add comment','Yes');
                else
                    choice = questdlg(['Is this multiplet of ' labb ' really a ' txt ' ?'], 'Validation of 1D multiplet', 'Yes','No','No, add comment','Yes');
                end
                opt.validate(1,3)=opt.validate(1,3)-1;
            end
        end
    end
    
    %set boudaries of axis
    margin1=0.2;fp1=power(10,-1);
    if abs(min_chem_shift_spectrum-max_chem_shift_spectrum)>10
        margin1=1;fp1=power(10,0);
    end
    if abs(min_chem_shift_spectrum-max_chem_shift_spectrum)>100
        margin1=10;fp1=power(10,1);
    end
    min_chem_shift_spectrum=min_chem_shift_spectrum -margin1;min_chem_shift_spectrum=round(min_chem_shift_spectrum/fp1)*fp1;
    max_chem_shift_spectrum=max_chem_shift_spectrum +margin1;max_chem_shift_spectrum=round(max_chem_shift_spectrum/fp1)*fp1;
    xlim([min_chem_shift_spectrum max_chem_shift_spectrum])
end
drawnow
end
