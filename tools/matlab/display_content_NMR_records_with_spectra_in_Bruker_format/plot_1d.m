function this_return=plot_1d(spectrum,hold_me,scales_in,sizes_in,fig_number)
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
end
if nargin<3
    scales_in=[min(spectrum.scale2) max(spectrum.scale2) ];
end
if nargin<2
    hold_me=0;
end
%% could include here interpolation...
    int_scale2=spectrum.scale2;
    data=spectrum.spectrum;
    


hhf=figure(fig_number);
if ~hold_me
    clf;
end
this_return=subplot(1,1,1);

hold on
box on


set(gca,'XDir','reverse');


% max(max(data));
% min(min(data));
%
% size(int_scale2)
% size(int_scale1)
% size(data)
plot(int_scale2,data)

drawnow


%axis([scales_in]);
%text(scales_in(1,2),scales_in(1,4),['F1/F2 ' num2str(round(factor1)) '/' num2str(round(factor2)) ])

end



