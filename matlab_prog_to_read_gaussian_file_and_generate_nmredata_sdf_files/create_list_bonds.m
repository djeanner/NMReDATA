function bond_list=create_list_bonds(atom_number,xyz,atoms_no_exchange,dist_max_xh,dist_max_cc,dist_max_cn,dist_max_co)
%% make list of bonds
chem_el='H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHzTlPbBiPoAtRnFrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDsRg';
%el=1;disp(['Element number ' num2str(el) ':' chem_el(el*2-1)  chem_el(el*2)]);% demo use chem_el


%% set distances for C-X bonds
disp_3d=1;%% to display bouds formed

dist_c_other=ones(3,120)*1000*0;
if nargin<7%for c-o
    dist_c_other(1,8)=1.8;
    dist_c_other(2,8)=1.3;%limit single/double
    %  dist_c_other(2,8)=1.3;
else
    dist_c_other(1,8)=dist_max_co;
end
if nargin<6 %for c-n
    dist_c_other(1,7)=sqrt(3);
else
    dist_c_other(1,7)=dist_max_cn;
end
if nargin<5 %for c-c
    dist_c_other(1,6)=1.8;
    dist_c_other(2,6)=1.44;
    %  dist_c_other(2,6)=1.5;
    dist_c_other(3,6)=1.25;
else
    dist_c_other(1,6)=dist_max_cc;
end
%% set distances for X-H bonds
if nargin<4
    dist_max_xh=sqrt(1.8); %maximal distance of CH bond
end

if nargin<3
    atoms_no_exchange=[6 7 8]; %with only 6, assume in dmso, with protons on N and O not exchanging...
    % atoms_no_exchange=[6]; %with only 6, assume in cdcl3, with protons on N and O exchanging...
end
if nargin<2
    xyz=[0 0 0 ; 0.5 0.5 0.5 ;-0.1 -0.1 -0.1]
end
if nargin<1
    atom_number=[6 ;6;1]
end
bond_list=zeros(size(xyz,1),size(xyz,1));

if disp_3d
    figure(888);clf;hold on;axis('equal');%pbaspect([1 1 1])
end
for looop_of_1=1:size(atom_number,1)%for each C ...
    text_tabel='';
    el=atom_number(looop_of_1,1);
    if el~=6
        text_tabel=[chem_el(el*2-1)  chem_el(el*2)];
    end
    text_tabel=[text_tabel '(' num2str(looop_of_1) ')'];
    text(xyz(looop_of_1,1) ,xyz(looop_of_1,2) ,xyz(looop_of_1,3) ,text_tabel,'HorizontalAlignment','center','VerticalAlignment','middle')
    %el=1;disp(['Element number ' num2str(el) ':' chem_el(el*2-1)  chem_el(el*2)]);% demo use chem_el
    
    if atom_number(looop_of_1,1)==6%if carbon
        
        for looop_of_other=1:size(atom_number,1)
            if looop_of_1~=looop_of_other
                for loop_over_het=atoms_no_exchange
                    
                    if atom_number(looop_of_other,1)==loop_over_het
                        tmp=xyz(looop_of_1,:)-xyz(looop_of_other,:);
                        dist=sqrt(sum(sum(tmp.*tmp)));
                        %% if single bound
                        widht_line=0;
                        if dist<dist_c_other(1,loop_over_het)
                            bond_list(looop_of_1,looop_of_other)=1;
                            bond_list(looop_of_other,looop_of_1)=1;
                            widht_line=1;line_type='k-';%for matlab drawing
                        end
                        %% if double bound
                        if dist<dist_c_other(2,loop_over_het)
                            bond_list(looop_of_1,looop_of_other)=2;
                            bond_list(looop_of_other,looop_of_1)=2;
                            widht_line=3;line_type='k-';%for matlab drawing
                            
                        end
                        %% if triple bound
                        if dist<dist_c_other(3,loop_over_het)
                            bond_list(looop_of_1,looop_of_other)=3;
                            bond_list(looop_of_other,looop_of_1)=3;
                            widht_line=3;line_type='k:';%for matlab drawing
                            
                        end
                        if disp_3d
                            if widht_line
                                plot3([xyz(looop_of_1,1) xyz(looop_of_other,1)],[xyz(looop_of_1,2) xyz(looop_of_other,2)],[xyz(looop_of_1,3) xyz(looop_of_other,3)],line_type,'LineWidth',widht_line)%3JHH
                            end
                        end
                    end
                end
            end
        end
    end
    if atom_number(looop_of_1,1)==1%if hydrogen
        for looop_of_other=1:size(atom_number,1)
            if looop_of_1~=looop_of_other
                if atom_number(looop_of_other,1)~=1,% no H-H
                    tmp=xyz(looop_of_1,:)-xyz(looop_of_other,:);
                    dist=sqrt(sum(sum(tmp.*tmp)));
                    if dist<dist_max_xh
                        bond_list(looop_of_1,looop_of_other)=1;
                        bond_list(looop_of_other,looop_of_1)=1;
                        if disp_3d
                            plot3([xyz(looop_of_1,1) xyz(looop_of_other,1)],[xyz(looop_of_1,2) xyz(looop_of_other,2)],[xyz(looop_of_1,3) xyz(looop_of_other,3)],'k:')%3JHH
                        end
                    end
                end
            end
        end
    end
end

end
