function this_return=plot_2d_interp(spectrum,hold_me,scales_in,sizes_in,fig_number)
%%

pt_per_mm=5;%resolution in pt per mm



if nargin<5
    fig_number=round(rand(1)*100000);
end

% plot region fast giving plot size
%eg plot_2d_interp(data,1,[],[])


% switch nargin
%     case 3
%         result = a + b;
%     case 1
%         result = a + a;
%     otherwise
%         result = 0;
% end
if nargin<4
    sizes_in=[0 0 15 10];
else
    if size(sizes_in,2)==0,
       
            sizes_in=[0 0 15 10];

    end
end
if nargin<3
    scales_in=[min(spectrum.scale2) max(spectrum.scale2) min(spectrum.scale1) max(spectrum.scale1)];
else
    if size(scales_in,2)==0
    scales_in=[min(spectrum.scale2) max(spectrum.scale2) min(spectrum.scale1) max(spectrum.scale1)];
    end
end
if nargin<2
    hold_me=0;
end
interpolate=1;

if interpolate
    %     pt_per_mm=5;%resolution in pt per mm
    %     nbpt2=pt_per_mm*10*abs(sizes_in(1,2)-sizes_in(1,1));
    %     nbpt1=pt_per_mm*10*abs(sizes_in(1,3)-sizes_in(1,4));
    %below keep all the same density of points... as first step.. then will
    %use interpft
    nbpt2=size(spectrum.spectrum,2)*abs(scales_in(1,2)-scales_in(1,1))/abs(min(spectrum.scale2)-max(spectrum.scale2));
    nbpt1=size(spectrum.spectrum,1)*abs(scales_in(1,3)-scales_in(1,4))/abs(min(spectrum.scale1)-max(spectrum.scale1));
    
    
    inc2=(abs(scales_in(1,2)-scales_in(1,1))    /(nbpt2-1) ) ;
    inc1=(abs(scales_in(1,3)-scales_in(1,4))    /(nbpt1-1) ) ;
    int_scale2=[scales_in(1,1):inc2:scales_in(1,2)];
    int_scale1=[scales_in(1,3):inc1:scales_in(1,4)];
    [mgr2, mgr1]=meshgrid(int_scale2,int_scale1);
    
    %     after=size(mgr2)
    %     before=size(data_set.spectrum)
    
    data=interp2(spectrum.scale2,spectrum.scale1,spectrum.spectrum,mgr2, mgr1,'nearest')   ;
    
    %now reduce number of points
%     sizes_in(1,3)
%     sizes_in(1,4)
%     pt_per_mm
    nbpt2=round(pt_per_mm*10*abs(sizes_in(1,3)));
    nbpt1=round(pt_per_mm*10*abs(sizes_in(1,4)));
    
    %round the factor
    %     factor2=round(size(data,2)/nbpt2)
    %     factor1=round(size(data,1)/nbpt1)
    %     nbpt1=round(size(data,1)/factor1)
    %     nbpt2=round(size(data,2)/factor2)
    factor2=(size(data,2)/nbpt2);
    factor1=(size(data,1)/nbpt1);
    if factor2<1,
        nbpt2=size(data,2);
        data2=data;
    else
        data2=zeros(size(data,1),nbpt2);
        for loo=1:nbpt2,
            from=round((loo-1)*factor2+1);
            to=round((loo-1)*factor2+factor2);
            %         if to>size(data,2),
            %             to=size(data,2)
            %             loo
            %             nbpt1
            %         end
            tmp=   data(:,[from:to]);
            tmp2=max(tmp,[],2);
            %             loo
            %             from
            %             to
            %             factor2
            %             factor1
            data2(:,loo)=tmp2;
        end
    end
    if factor1<1,
        nbpt1=size(data,1);
        data=data2;
    else
        
        data=zeros(nbpt1,nbpt2);
        for loo=1:nbpt1,
            % (loo-1)*factor1+[1:factor1];
            from=round((loo-1)*factor1+1);
            to=round((loo-1)*factor1+factor1);
            tmp=   data2([from:to],:);
            tmp2=max(tmp,[],1);
            data(loo,:)=tmp2;
        end
    end
    inc2=(abs(scales_in(1,2)-scales_in(1,1))    /(nbpt2-1) ) ;
    inc1=(abs(scales_in(1,3)-scales_in(1,4))    /(nbpt1-1) ) ;
    int_scale2=[scales_in(1,1):inc2:scales_in(1,2)];
    int_scale1=[scales_in(1,3):inc1:scales_in(1,4)];
    
    %      toc
    %      data=interpft(data,nbpt1,1);
    %      data(isnan(data))=0;
    %
    %     toc
    %      data=interpft(data,nbpt2,2);
    %      data(isnan(data))=0;
    
    
    
else
    int_scale2=spectrum.scale2;
    int_scale1=spectrum.scale1;
    data=spectrum.spectrum;
    
end

cont_level_list= spectrum.cont_level_list;
hhf=figure(fig_number);
if ~hold_me,
    clf;
end
this_return=subplot(1,1,1);

hold on
box on


set(gca,'XDir','reverse');
set(gca,'YDir','reverse');

max_list_cont=size(cont_level_list,2);
my_col_map1=colormap(autumn(max_list_cont+1));%pos contours
my_col_map2=colormap(summer(max_list_cont+1));%neg contours


%black and yellow
my_col_map1=colormap(gray(max_list_cont+1));%pos contours
my_col_map2=colormap(pink(max_list_cont+1));%neg contours

%changed to blue/red....
my_col_map2=[ones(size(my_col_map1,1),1)*0001,ones(size(my_col_map1,1),1)*0000,ones(size(my_col_map1,1),1)*00000];
my_col_map1=[ones(size(my_col_map1,1),1)*0000,ones(size(my_col_map1,1),1)*0000,ones(size(my_col_map1,1),1)*00001];



% max(max(data));
% min(min(data));
%
% size(int_scale2)
% size(int_scale1)
% size(data)
for list_contour=1:max_list_cont,
    now= cont_level_list(list_contour);
    contour(int_scale2,int_scale1,data,[   cont_level_list(list_contour)  cont_level_list(list_contour)],'Linecolor',my_col_map1(max_list_cont-list_contour+1,:));
    contour(int_scale2,int_scale1,data,[  -cont_level_list(list_contour) -cont_level_list(list_contour)],'Linecolor',my_col_map2(max_list_cont-list_contour+1,:));
    drawnow;
    
end
drawnow

drawnow;
if ~interpolate,
    axis([scales_in]);
end

%axis([scales_in]);
%text(scales_in(1,2),scales_in(1,4),['F1/F2 ' num2str(round(factor1)) '/' num2str(round(factor2)) ])
text_cont=' ';
if round(factor1)>1
    text_cont=[text_cont 'F1 1/' num2str(round(factor1)) ' '  ];
end
if round(factor2)>1
    text_cont=[text_cont 'F2 1/' num2str(round(factor2)) ' '  ];
end
    text_cont=[text_cont ' ' num2str(pt_per_mm) ' pt/mm'  ];

text(scales_in(1,2),scales_in(1,4),text_cont,'FontSize',8,'VerticalAlignment','bottom')
end



