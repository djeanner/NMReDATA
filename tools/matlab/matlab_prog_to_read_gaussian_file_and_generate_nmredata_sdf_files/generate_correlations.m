function [cosy, hsqc, hmbc]=generate_correlations(system,atoms_no_exchange,disp_3d)
atom_number=system.atom_number;
bond_list=system.bond_list;

if nargin<3
    disp_3d=1;
end
if isfield(system,'xyz')
xyz=system.xyz;
else
        disp_3d=0;

end
%function [ cosy, hsqc, hmbc]=read_dft_nmr_data_fn(dft_folder,min_J)
%only for C-H (no OH !)
disp_cosy=1;
disp_hmbc=1;



if nargin<2
    atoms_no_exchange=[6 7 8]; %with only 6, assume in dmso, with protons on N and O not exchanging...
    atoms_no_exchange=[6]; %with only 6, assume in cdcl3, with protons on N and O exchanging...
end

hsqc=zeros(size(atom_number,1),size(atom_number,1));
cosy=hsqc*0;

number_of_h=zeros(1,size(atom_number,1));
for looop_of_cno=1:size(atom_number,1)%for each C ...
    for loop_over_het=atoms_no_exchange
        if atom_number(looop_of_cno,1)==loop_over_het%if C/O/N
            for looop_of_h=1:size(atom_number,1)
                if atom_number(looop_of_h,1)==1% if H...
                    % tmp=xyz(looop_of_c,:)-xyz(looop_of_h,:);
                    % dist=sum(sum(tmp.*tmp));
                    if bond_list(looop_of_cno,looop_of_h)
                        if disp_3d
                            plot3([xyz(looop_of_cno,1) xyz(looop_of_h,1)],[xyz(looop_of_cno,2) xyz(looop_of_h,2)],[xyz(looop_of_cno,3) xyz(looop_of_h,3)],'b-')%CH
                        end
                       % if atom_number(loop_over_het,1)==6
                        if loop_over_het==6
                            hsqc(looop_of_h,looop_of_cno)=1;
                            hsqc(looop_of_cno,looop_of_h)=1;
                        end
                        %% cosy 2J
                        for loop_over_h2=1:size(atom_number,1)
                            if (atom_number(loop_over_h2,1)==1) && (loop_over_h2~=looop_of_h) && hsqc(looop_of_cno,loop_over_h2)
                                cosy(loop_over_h2,looop_of_h)=2;
                                cosy(looop_of_h,loop_over_h2)=2;
                                if disp_3d && disp_cosy
                                    plot3([xyz(looop_of_h,1) xyz(loop_over_h2,1)],[xyz(looop_of_h,2) xyz(loop_over_h2,2)],[xyz(looop_of_h,3) xyz(loop_over_h2,3)],'m:');%2JHH
                                    text([xyz(looop_of_h,1) *0.5+0.5* xyz(loop_over_h2,1)],[xyz(looop_of_h,2) *0.5+0.5* xyz(loop_over_h2,2)],[xyz(looop_of_h,3) *0.5+0.5* xyz(loop_over_h2,3)],num2str(system.J(loop_over_h2,looop_of_h)));%2JHH
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% list carbons connected to carbons (for list of HMBC correlations)
hmbc=hsqc*0;

for looop_of_cno=1:size(atom_number,1)%for each C ...
    for loop_over_het=atoms_no_exchange
        if atom_number(looop_of_cno,1)==loop_over_het%if carbon
            list_c{looop_of_cno}=[];
            count=1;
            for looop_of_c2=1:size(atom_number,1)
                for loop_over_het2=atoms_no_exchange
                    if (atom_number(looop_of_c2,1)==loop_over_het2) && (looop_of_c2~=looop_of_cno)% if C...
                        % tmp=xyz(looop_of_c,:)-xyz(looop_of_c2,:);
                        % dist=sum(sum(tmp.*tmp));
                        % if dist<dist_max_cc
                        if bond_list(looop_of_cno,looop_of_c2)
                            list_c{looop_of_cno}=[list_c{looop_of_cno} looop_of_c2];
                            if disp_3d
                                plot3([xyz(looop_of_cno,1) xyz(looop_of_c2,1)],[xyz(looop_of_cno,2) xyz(looop_of_c2,2)],[xyz(looop_of_cno,3) xyz(looop_of_c2,3)],'k-')% CC
                            end
                            for loop_over_h1=1:size(atom_number,1)%
                                if atom_number(loop_over_h1,1)==1
                                    if bond_list(loop_over_h1,looop_of_cno)
                                        if atom_number(looop_of_c2,1)==6
                                            
                                            hmbc(loop_over_h1,looop_of_c2)=2;
                                            hmbc(looop_of_c2,loop_over_h1)=2;
                                        end
                                        if disp_3d && disp_hmbc
                                            plot3([xyz(loop_over_h1,1) xyz(looop_of_c2,1)],[xyz(loop_over_h1,2) xyz(looop_of_c2,2)],[xyz(loop_over_h1,3) xyz(looop_of_c2,3)],'c:')% 2JCH
                                            text([xyz(loop_over_h1,1) *0.5+0.5*xyz(looop_of_c2,1)],[xyz(loop_over_h1,2)*0.5+0.5* xyz(looop_of_c2,2)],[xyz(loop_over_h1,3) *0.5+0.5*xyz(looop_of_c2,3)],num2str(system.J(loop_over_h1,looop_of_c2)))% 2JCH
                                        end
                                        for loop_over_h2=1:size(atom_number,1)%
                                            if atom_number(loop_over_h2,1)==1
                                                if bond_list(loop_over_h2,looop_of_c2)
                                                    cosy(loop_over_h2,loop_over_h1)=3;
                                                    cosy(loop_over_h1,loop_over_h2)=3;
                                                    if disp_3d && disp_cosy
                                                        plot3([xyz(loop_over_h1,1) xyz(loop_over_h2,1)],[xyz(loop_over_h1,2) xyz(loop_over_h2,2)],[xyz(loop_over_h1,3) xyz(loop_over_h2,3)],'m--')%3JHH
                                                        text([xyz(loop_over_h1,1)*0.5+0.5* xyz(loop_over_h2,1)],[xyz(loop_over_h1,2) *0.5+0.5*xyz(loop_over_h2,2)],[xyz(loop_over_h1,3)*0.5+0.5* xyz(loop_over_h2,3)],num2str(system.J(loop_over_h2,loop_over_h1)))%3JHH
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
% similar as above...
for looop_of_cno=1:size(atom_number,1)%for each C ...
    for loop_over_het=atoms_no_exchange
        
        if atom_number(looop_of_cno,1)==6%if carbon
            list_c{looop_of_cno}=[];
            count=1;
            for looop_of_c2=1:size(atom_number,1)
                if (atom_number(looop_of_c2,1)==6) && (looop_of_c2~=looop_of_cno)% if C...
                    % tmp=xyz(looop_of_c,:)-xyz(looop_of_c2,:);
                    %   dist=sum(sum(tmp.*tmp));
                    %  if dist<dist_max_cc
                    if bond_list(looop_of_cno,looop_of_c2)
                        list_c{looop_of_cno}=[list_c{looop_of_cno} looop_of_c2];
                        if disp_3d
                            plot3([xyz(looop_of_cno,1) xyz(looop_of_c2,1)],[xyz(looop_of_cno,2) xyz(looop_of_c2,2)],[xyz(looop_of_cno,3) xyz(looop_of_c2,3)],'k-')% CC
                        end
                        for loop_over_h1=1:size(atom_number,1)%
                            if atom_number(loop_over_h1,1)==1
                                if bond_list(loop_over_h1,looop_of_cno)
                                    %%for 3JCH hmbc
                                    if size(list_c,2)>=looop_of_c2
                                        if size(list_c{looop_of_c2},2)>0
                                            for loop_over_next_c=list_c{looop_of_c2}
                                                if loop_over_next_c~=looop_of_cno
                                                    if atom_number(loop_over_next_c,1)==6
                                                        hmbc(loop_over_h1,loop_over_next_c)=3;
                                                        hmbc(loop_over_next_c,loop_over_h1)=3;
                                                    end
                                                    if disp_3d && disp_hmbc
                                                        plot3([xyz(loop_over_h1,1) xyz(loop_over_next_c,1)],[xyz(loop_over_h1,2) xyz(loop_over_next_c,2)],[xyz(loop_over_h1,3) xyz(loop_over_next_c,3)],'c--')%3JCH
                                                        text([xyz(loop_over_h1,1)*0.5+0.5* xyz(loop_over_next_c,1)],[xyz(loop_over_h1,2)*0.5+0.5* xyz(loop_over_next_c,2)],[xyz(loop_over_h1,3) *0.5+0.5*xyz(loop_over_next_c,3)],num2str(system.J(loop_over_h1,loop_over_next_c)))%3JCH
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end
