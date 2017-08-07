function draw_coupling_plot_new(moulin_dft,red_system ,min_coupling_for_display)
if nargin<3,
    min_coupling_for_display=0.5;
end

nucleus_considered1=1;%
nucleus_considered2=1;%
%%%type_of_bound=-1;%-1 for all number of bonds, 1 for 1, 2, for 2, 3, for 3, 23 for two and threee 2345 for two and more.
% not working...

list1=find(red_system.atom_number==nucleus_considered1)';
list2=find(red_system.atom_number==nucleus_considered2)';
table_cs=red_system.cs(:,1)';

mi=min(min(table_cs(1,list1)));
ma=max(max(table_cs(1,list1)));

lin_spaced=0;
if lin_spaced
    table_cs2=(mi:(ma-mi)/(size(table_cs(1,list1),2)-1):ma);
    table_cs(1,list1)=table_cs2;
end

%
%
%                 xyz: [43?3 double]
%         atom_number: [43?1 int32]
%                   J: [43?43 double]
%                  cs: [43?1 double]
%           bond_list: [43?43 double]
%         label_atoms: {1?43 cell}
%              lookup: {1?43 cell}
%     nb_bond_between: [43?43 double]
%


figure(33);clf;hold on
y_pos=moulin_dft.from_y-moulin_dft.width_y*min_coupling_for_display/moulin_dft.max_cou;
plot([mi-0.1 ma+0.1],[y_pos y_pos],'r:')

for inc=list1
    % for inc=1:size(coupling_table,1),
    
    % plot support for dots
    ticks_pos=0;
    plot(table_cs(1,inc), moulin_dft.from_y-moulin_dft.width_y*ticks_pos/moulin_dft.max_cou,'k.');%origin of line
    
    for ticks_pos=moulin_dft.ticks_at:moulin_dft.ticks_at:moulin_dft.max_cou-0.1,
        plot(table_cs(1,inc), moulin_dft.from_y-moulin_dft.width_y*ticks_pos/moulin_dft.max_cou,'k+');% main ticksline
    end
    ticks_pos=moulin_dft.max_cou;
    plot(table_cs(1,inc), moulin_dft.from_y-moulin_dft.width_y*ticks_pos/moulin_dft.max_cou,'k^');%top of arrow line
    ticks_pos=moulin_dft.max_cou+0.2*moulin_dft.ticks_at;
    plot(table_cs(1,inc)*[1 1], [moulin_dft.from_y moulin_dft.from_y-moulin_dft.width_y],'k-');%main line
    
    % add label
    txt=[ '' red_system.label_atoms{inc}];
    this_refffj= text(table_cs(1,inc), moulin_dft.from_y-moulin_dft.width_y*ticks_pos/moulin_dft.max_cou,txt);%line
    set(this_refffj,'rotation',90);
    for loj=list2
        if loj~=inc
            if abs(red_system.J(inc,loj))>min_coupling_for_display
                if round(1+127*abs(red_system.J(inc,loj))/moulin_dft.max_cou)>0 &&  round(round(1+127*abs(red_system.J(inc,loj))/moulin_dft.max_cou))<=size(moulin_dft.my_color_map,1)
                    my_col=moulin_dft.my_color_map(round(1+127*abs(red_system.J(inc,loj))/moulin_dft.max_cou),:);
                else
                    my_col=[0 0 0];
                end
                loopgetpartner=loj;
                y_pos=moulin_dft.from_y-moulin_dft.width_y*abs(red_system.J(inc,loj))/moulin_dft.max_cou;
                ref_for_c=plot(table_cs(1,inc), y_pos,'o','color',my_col);%point
                ref_for_c2=plot([table_cs(1,inc) table_cs(1,loj) ], [y_pos y_pos] ,'-');%point
                ref_for_c3=text([table_cs(1,inc)*0.5+0.5* table_cs(1,loopgetpartner) ], y_pos,num2str(red_system.J(inc, loopgetpartner),'%.1f'),'HorizontalAlignment','center');%line
                set(ref_for_c2,'Color',my_col);
                set(ref_for_c3,'Color',my_col);
                
            end
        end
    end
end
set(gca,'Xdir','reverse')
ylim([0 26])
