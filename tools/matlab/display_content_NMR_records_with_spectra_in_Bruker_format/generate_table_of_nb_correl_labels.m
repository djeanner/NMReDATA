function  [nb_bond_between_atoms_including_implicit_H,ref_from_label_number_to_nb_bond_table]=generate_table_of_nb_correl_labels(object)
verbose=1;
nb_bond_between_atoms_including_implicit_H=object.structure.nb_bond_between;% create_nb_bonds_between_atoms(structure,opt);
ref_from_label_number_to_nb_bond_table=zeros(1,size(object.atom_num,2));

%first run for explicit atoms
for lo_over_labels=1:size(object.label_signal,2)
    last_tmp=0;
    
    if verbose
        % disp(['For label : ' object.label_signal{1,lo_over_labels}]);
    end
    counter=0;
    
    for lo_over_atoms=1:size(object.atom_num,2)
        tmp=object.atom_num{lo_over_labels,lo_over_atoms};
        if lo_over_atoms==2
            %  disp(['multiple atoms to this label' tmp])
        end
        if size(tmp,2)>0
            counter=counter+1;
            if verbose
                disp(['Label ' num2str(lo_over_labels) ': "' object.label_signal{1,lo_over_labels} '" atom(' num2str(lo_over_atoms) '):' num2str( object.atom_num{lo_over_labels,lo_over_atoms})])
            end
            last_tmp=tmp;
        end
    end
    %  disp(['counter:' num2str(counter)])
    if counter==1
        if last_tmp>0
            ref_from_label_number_to_nb_bond_table(lo_over_labels)=last_tmp;
        end
        
    end
end

%second run to work on nb_bond_between_H table (add implicit H...)
for lo_over_labels=1:size(object.label_signal,2)
    last_tmp=0;
    
    if verbose
        % disp(['For label : ' object.label_signal{1,lo_over_labels}]);
    end
    counter=0;
    
    for lo_over_atoms=1:size(object.atom_num,2)
        tmp=object.atom_num{lo_over_labels,lo_over_atoms};
        if lo_over_atoms==2
            %  disp(['multiple atoms to this label' tmp])
        end
        if size(tmp,2)>0
            counter=counter+1;
            
            last_tmp=tmp;
        end
    end
    %  disp(['counter:' num2str(counter)])
    if counter==1
        if last_tmp<0
            if verbose
                disp(['Implicit H for atom ' num2str(-last_tmp) ' for label ' object.label_signal{1,lo_over_labels} '. Inserts a line number ' num2str(size(nb_bond_between_atoms_including_implicit_H,2)+1) ' in nb_bond_between_H table'])
                point_tmp=(-last_tmp);
                extract1= nb_bond_between_atoms_including_implicit_H(:,point_tmp,:);
                extract1(2:end,:,:)=extract1(1:end-1,:,:);%shift one step because all atoms are one bound further away
                extract1(1,size(extract1,2),:)=0*extract1(1,size(extract1,2),:);%set line to zero
                extract1(1,size(extract1,2),point_tmp)=1;%set directly attached atom one bond away
                nb_bond_between_atoms_including_implicit_H(:,size(nb_bond_between_atoms_including_implicit_H,3)+1,:)=extract1(:,size(extract1,2),:);
                extract2= nb_bond_between_atoms_including_implicit_H(:,:,point_tmp);
                extract2(2:end,:)=extract2(1:end-1,:);
                extract2(1,:)=0*extract2(1,:);
                extract2(1,point_tmp)=1;
                nb_bond_between_atoms_including_implicit_H(:,:,size(nb_bond_between_atoms_including_implicit_H,3)+1)=extract2;
                %store pointer:
                ref_from_label_number_to_nb_bond_table(lo_over_labels)=size(nb_bond_between_atoms_including_implicit_H,2);
            end
        end
    end
end
if verbose
    disp(['end ' ])
end

end