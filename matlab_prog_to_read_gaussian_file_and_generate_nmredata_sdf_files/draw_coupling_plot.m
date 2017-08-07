function draw_coupling_plot(moulin_dft,table_1h_cs,table_coup_for_return,table_Jxy,label_atoms,ref_to_atom)


figure(3);clf;hold on

for inc=1:size(table_coup_for_return,1)
    % for inc=1:size(coupling_table,1),
    
    
    ticks_pos=0;
    plot(table_1h_cs(1,inc), moulin_dft.from_y-moulin_dft.width_y*ticks_pos/moulin_dft.max_cou,'k.');%line
    
    for ticks_pos=moulin_dft.ticks_at:moulin_dft.ticks_at:moulin_dft.max_cou-0.1,
        plot(table_1h_cs(1,inc), moulin_dft.from_y-moulin_dft.width_y*ticks_pos/moulin_dft.max_cou,'k+');%line
    end
    ticks_pos=moulin_dft.max_cou;
    plot(table_1h_cs(1,inc), moulin_dft.from_y-moulin_dft.width_y*ticks_pos/moulin_dft.max_cou,'k^');%line
    ticks_pos=moulin_dft.max_cou+0.2*moulin_dft.ticks_at;
    %
    txt='';
    if size(ref_to_atom,1)>=inc
    txt=[ '' label_atoms(ref_to_atom(inc))];
    end
    
    this_refffj= text(table_1h_cs(1,inc), moulin_dft.from_y-moulin_dft.width_y*ticks_pos/moulin_dft.max_cou,txt);%line
    set(this_refffj,'rotation',90);
    
    
    plot(table_1h_cs(1,inc)*[1 1], [moulin_dft.from_y moulin_dft.from_y-moulin_dft.width_y],'k-');%line
    for loj=1:size(table_coup_for_return,2)
        if table_coup_for_return(inc,loj)==0  
            break;
        else
            %%   if coupling_table(inc,loj)<moulin_dft.max_cou,%othersise error in colormap expno
            
            ref_for_c=plot(table_1h_cs(1,inc), moulin_dft.from_y-moulin_dft.width_y*table_coup_for_return(inc,loj)/moulin_dft.max_cou,'o');%point
            for loopgetpartner=inc+1:size(table_Jxy ,1)
                if abs(table_coup_for_return(inc,loj)-abs(table_Jxy(inc, loopgetpartner)))<0.00001,% table_Jxy coup max be negative second abs
                    if table_Jxy(inc, loopgetpartner)>0
                        ref_for_c2=plot([table_1h_cs(1,inc) table_1h_cs(1,loopgetpartner) ], [1 1]*moulin_dft.from_y-moulin_dft.width_y*table_coup_for_return(inc,loj)/moulin_dft.max_cou,':');%line
                    else
                        ref_for_c2=plot([table_1h_cs(1,inc) table_1h_cs(1,loopgetpartner) ], [1 1]*moulin_dft.from_y-moulin_dft.width_y*table_coup_for_return(inc,loj)/moulin_dft.max_cou,'--');%line
                    end
                    ref_for_c3=text([table_1h_cs(1,inc)*0.5+0.5* table_1h_cs(1,loopgetpartner) ], moulin_dft.from_y-moulin_dft.width_y*table_coup_for_return(inc,loj)/moulin_dft.max_cou,num2str(table_Jxy(inc, loopgetpartner),'%.1f'));%line
                                if round(1+127*table_coup_for_return(inc,loj))<=size(moulin_dft.my_color_map)
                    set(ref_for_c3,'Color',moulin_dft.my_color_map(round(1+127*table_coup_for_return(inc,loj)/moulin_dft.max_cou),:));
                    set(ref_for_c2,'Color',moulin_dft.my_color_map(round(1+127*table_coup_for_return(inc,loj)/moulin_dft.max_cou),:));
                                end
                end
            end
            if round(1+127*table_coup_for_return(inc,loj))<=size(moulin_dft.my_color_map)
            set(ref_for_c,'Color',moulin_dft.my_color_map(round(1+127*table_coup_for_return(inc,loj)/moulin_dft.max_cou),:));
            end
            %   set(ref_for_c,'Color',moulin_dft.my_color_map(round(1+127*table_coup_for_return(inc,loj)/moulin_dft.max_cou),:));
            
            %%  end
        end
    end
end
set(gca,'Xdir','reverse')
ylim([0 26])