function export_to_acd(nmredata)

%% opening file
fileID = fopen('demo_exported_ACD.sdf','w');

%% write the mol block
mol_block=nmredata{1}.structure.mol_block;
for loop=1:size(mol_block,2)
    fprintf(fileID,'%s\n',nmredata{1}.structure.mol_block{1,loop});
end

%% scan all spectra
for loop_over_elements=1:size(nmredata,2)
    cur_obj=nmredata{1,loop_over_elements};
    if isfield(cur_obj,'nb_dim')
        if cur_obj.nb_dim==1
            disp(['Found a 1D spectrum (object = ' num2str(loop_over_elements) ') : ' cur_obj.tag_name])
            disp(['list of fields for this object:' ])
            disp(['-------------------------------------------------' ])
            cur_obj
            disp(['-------------------------------------------------' ])
            
            if isfield(cur_obj,'chemical_shift') && isfield(cur_obj,'label') 
                %% write tag
                fprintf(fileID,'>  <DEMO_TAG_%d>\n',loop_over_elements);%start of tag
                fprintf(fileID,'%s : %f ppm\n',cur_obj.label{1},cur_obj.chemical_shift{1});%start of tag
                fprintf(fileID,'\n',loop_over_elements);%empty line as tag separator
            end
            
        end
        if cur_obj.nb_dim==2
            disp(['Found a 2D spectrum (object = ' num2str(loop_over_elements) ') : ' cur_obj.tag_name])
            disp(['list of fields for this object:' ])
            disp(['-------------------------------------------------' ])
            cur_obj
            disp(['-------------------------------------------------' ])
        end
    end
end


%% close file
fclose(fileID)


end