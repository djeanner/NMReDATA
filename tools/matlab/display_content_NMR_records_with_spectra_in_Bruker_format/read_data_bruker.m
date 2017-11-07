function data=read_data_bruker(path_short,acquno,procno,all_quadrants)
%% New 9.3.2016 read all d and cnst
%% corrected intensity scale... ncproc
%% reads xf2 data
if (nargin<4)
    all_quadrants=0;
end
data=struct;

k = strfind(path_short,'/');
data.path_base='/';
if size(k,2)>1
    data.path_base=path_short(1:k(end-1));
end


data.path=path_short;
data.acquno=acquno;
data.procno=procno;
data.full_data=[ num2str(acquno)  '/' num2str(procno) ];
list_slash=strfind(data.path,'/');
l2=list_slash(1,end);
l1=list_slash(1,end-1);

data.expname=data.path((l1+1):(l2-1));
data.expnamefullshort=[data.expname '/' num2str(acquno)  '/' num2str(procno)];

si1=1;
path=[path_short num2str(acquno) '/'];
if exist([path '/acqus'],'file')
    data.procno=procno;
    data.acquno=acquno;
    content=textread([path '/acqus'],'%s');
    data.rg         = str2num(char(content(strmatch('##$RG=',content')+1)));
    data.ns         = str2num(char(content(strmatch('##$NS=',content')+1)));
    data.td2         = str2num(char(content(strmatch('##$TD=',content')+1)));
    data.o2         = str2num(char(content(strmatch('##$O2=',content')+1)));
    data.sfo2         = str2num(char(content(strmatch('##$SFO2=',content')+1)));
    data.bf2         = str2num(char(content(strmatch('##$BF2=',content')+1)));
    data.o1         = str2num(char(content(strmatch('##$O1=',content')+1)));
    data.sfo1         = str2num(char(content(strmatch('##$SFO1=',content')+1)));
    data.bf1         = str2num(char(content(strmatch('##$BF1=',content')+1)));
    data.pulprog         = (char(content(strmatch('##$PULPROG=',content')+1)));
    data.o1p=data.o1/data.bf1;
    data.o2p=data.o2/data.bf2;
    data.cnst30         =str2num(char(content(strmatch('##$CNST=',content')+2+30)));
    data.d7         =str2num(char(content(strmatch('##$D=',content')+2+7)));
    data.d8         =str2num(char(content(strmatch('##$D=',content')+2+8)));
    data.d9         =str2num(char(content(strmatch('##$D=',content')+2+9)));
    data.d10         =str2num(char(content(strmatch('##$D=',content')+2+10)));
    data.d11         =str2num(char(content(strmatch('##$D=',content')+2+11)));
    data.d12         =str2num(char(content(strmatch('##$D=',content')+2+12)));
    data.d13         =str2num(char(content(strmatch('##$D=',content')+2+13)));
    data.d14         =str2num(char(content(strmatch('##$D=',content')+2+14)));
    data.d15         =str2num(char(content(strmatch('##$D=',content')+2+15)));
    data.d16         =str2num(char(content(strmatch('##$D=',content')+2+16)));
    data.d17         =str2num(char(content(strmatch('##$D=',content')+2+17)));
    data.d18         =str2num(char(content(strmatch('##$D=',content')+2+18)));
    data.d19         =str2num(char(content(strmatch('##$D=',content')+2+19)));
    data.d20         =str2num(char(content(strmatch('##$D=',content')+2+20)));
    data.d21         =str2num(char(content(strmatch('##$D=',content')+2+21)));
    data.d22         =str2num(char(content(strmatch('##$D=',content')+2+22)));
    data.d23         =str2num(char(content(strmatch('##$D=',content')+2+23)));
    data.d24         =str2num(char(content(strmatch('##$D=',content')+2+24)));
    data.d25         =str2num(char(content(strmatch('##$D=',content')+2+25)));
    data.d26         =str2num(char(content(strmatch('##$D=',content')+2+26)));
    data.d27         =str2num(char(content(strmatch('##$D=',content')+2+27)));
    data.d28         =str2num(char(content(strmatch('##$D=',content')+2+28)));
    data.d29         =str2num(char(content(strmatch('##$D=',content')+2+29)));
    data.d30         =str2num(char(content(strmatch('##$D=',content')+2+30)));
    data.d31         =str2num(char(content(strmatch('##$D=',content')+2+31)));
    data.cnst1          =str2num(char(content(strmatch('##$CNST=',content')+2+1)));
    data.cnst2          =str2num(char(content(strmatch('##$CNST=',content')+2+2)));
    data.cnst3          =str2num(char(content(strmatch('##$CNST=',content')+2+3)));
    data.cnst4          =str2num(char(content(strmatch('##$CNST=',content')+2+4)));
    data.cnst5          =str2num(char(content(strmatch('##$CNST=',content')+2+5)));
    data.cnst6          =str2num(char(content(strmatch('##$CNST=',content')+2+6)));
    data.cnst7          =str2num(char(content(strmatch('##$CNST=',content')+2+7)));
    data.cnst8          =str2num(char(content(strmatch('##$CNST=',content')+2+8)));
    data.cnst9          =str2num(char(content(strmatch('##$CNST=',content')+2+9)));
    data.cnst10          =str2num(char(content(strmatch('##$CNST=',content')+2+10)));
    data.cnst11          =str2num(char(content(strmatch('##$CNST=',content')+2+11)));
    data.cnst12          =str2num(char(content(strmatch('##$CNST=',content')+2+12)));
    data.cnst13          =str2num(char(content(strmatch('##$CNST=',content')+2+13)));
    data.cnst14          =str2num(char(content(strmatch('##$CNST=',content')+2+14)));
    data.cnst15          =str2num(char(content(strmatch('##$CNST=',content')+2+15)));
    data.cnst16          =str2num(char(content(strmatch('##$CNST=',content')+2+16)));
    data.cnst17          =str2num(char(content(strmatch('##$CNST=',content')+2+17)));
    data.cnst18          =str2num(char(content(strmatch('##$CNST=',content')+2+18)));
    data.cnst19          =str2num(char(content(strmatch('##$CNST=',content')+2+19)));
    data.cnst20          =str2num(char(content(strmatch('##$CNST=',content')+2+20)));
    data.cnst21          =str2num(char(content(strmatch('##$CNST=',content')+2+21)));
    data.cnst22          =str2num(char(content(strmatch('##$CNST=',content')+2+22)));
    data.cnst23          =str2num(char(content(strmatch('##$CNST=',content')+2+23)));
    data.cnst24          =str2num(char(content(strmatch('##$CNST=',content')+2+24)));
    data.cnst25          =str2num(char(content(strmatch('##$CNST=',content')+2+25)));
    data.cnst26          =str2num(char(content(strmatch('##$CNST=',content')+2+26)));
    data.cnst27          =str2num(char(content(strmatch('##$CNST=',content')+2+27)));
    data.cnst28          =str2num(char(content(strmatch('##$CNST=',content')+2+28)));
    data.cnst29          =str2num(char(content(strmatch('##$CNST=',content')+2+29)));
    data.cnst30          =str2num(char(content(strmatch('##$CNST=',content')+2+30)));
    data.cnst31          =str2num(char(content(strmatch('##$CNST=',content')+2+31)));
    data.sw2= str2num(char(content(strmatch('##$SW=',content')+1)));
    data.sw2h= str2num(char(content(strmatch('##$SW_h=',content')+1)));
        
else
    data.file_not_found=[path '/acqus'];
    disp(['file not found : ' data.file_not_found]);
end
if exist([path '/acqu2s'])
    % read accq2s
    content=textread([path '/acqu2s'],'%s');
    data.fnmode         = str2num(char(content(strmatch('##$FnMODE=',content')+1)));
    data.td1         = str2num(char(content(strmatch('##$TD=',content')+1)));
    data.sw1         = str2num(char(content(strmatch('##$SW=',content')+1)));
    data.sw1h         = str2num(char(content(strmatch('##$SW_h=',content')+1)));
    
    % read procs
end

if exist([path '/pdata/' num2str(procno) '/procs'],'file')
    content=textread([path '/pdata/' num2str(procno) '/procs'],'%s');
    data.si2         = str2num(char(content(strmatch('##$SI=',content')+1)));
    data.wdw         = str2num(char(content(strmatch('##$WDW=',content')+1)));
    data.ssb         = str2num(char(content(strmatch('##$SSB=',content')+1)));
    data.phc0         = str2num(char(content(strmatch('##$PHC0=',content')+1)));
    data.phc1         = str2num(char(content(strmatch('##$PHC1=',content')+1)));
    data.xdim         = str2num(char(content(strmatch('##$XDIM=',content')+1)));
    data.tdeff2         = str2num(char(content(strmatch('##$TDeff=',content')+1)));
        data.sf2         = str2num(char(content(strmatch('##$SF=',content')+1)));

    data.scale2=[data.o1p+data.sw2/2:-data.sw2/data.si2:(data.o1p-data.sw2/2)+data.sw2/data.si2];
    
    si2=data.si2;
    nc_proc_tmp         = str2num(char(content(strmatch('##$NC_proc=',content')+1)));
    %NC_diff=-14-NC_proc; was this in previous versions...
    nc_diff=-nc_proc_tmp;
    nc=power(2,nc_diff);
    data.nc_proc=nc;
    data.NC_proc=nc_proc_tmp;%orig value
end

     
if exist([path '/pdata/' num2str(procno) '/proc2s'],'file')
    
    content=textread([path '/pdata/' num2str(procno) '/proc2s'],'%s');
    data.xdim1         = str2num(char(content(strmatch('##$XDIM=',content')+1)));
    data.si1         = str2num(char(content(strmatch('##$SI=',content')+1)));
    data.wdw1         = str2num(char(content(strmatch('##$WDW=',content')+1)));
    data.ssb1         = str2num(char(content(strmatch('##$SSB=',content')+1)));
    data.phc01         = str2num(char(content(strmatch('##$PHC0=',content')+1)));
    data.phc11         = str2num(char(content(strmatch('##$PHC1=',content')+1)));
    data.ftsize         = str2num(char(content(strmatch('##$FTSIZE=',content')+1)));
    data.tdeff1         = str2num(char(content(strmatch('##$TDeff=',content')+1)));
    data.swproc1         = str2num(char(content(strmatch('##$SW_p=',content')+1)));
    data.sfproc1        = str2num(char(content(strmatch('##$SF=',content')+1)));
    data.offset1         = str2num(char(content(strmatch('##$OFFSET=',content')+1)));
    data.sf1         = str2num(char(content(strmatch('##$SF=',content')+1)));
    data.hzppt1=data.swproc1/data.si1;
    
    data.swprocppm1=data.hzppt1*data.si1/data.sfproc1;
    
    %     NC_proc         = str2num(char(content(strmatch('##$NC_proc=',content')+1)));
    %
    %     NC_diff=-NC_proc;
    %     nc=power(2,NC_diff);
    %     data.nc_proc=nc;
    sw=data.swprocppm1;
    opk=data.offset1-sw/2;
    data.o2p_proc_1=opk;
    data.scale1=[opk+sw/2:-sw/data.si1:(opk-sw/2)+sw/data.si1];
    si1=data.si1;
end

%determines center in f1 (o1p or o2p)
if isfield(data,'ftsize')
    if data.ftsize>0%only if FT data in f1
        if data.sf1==data.sf2
            disp(['assuming homonuclear experiment set middle to o1p '])
            middle_of_spectrum_in_f1=data.o1p;
        else
            disp(['assuming heteronuclear experiment set middle to o2p '])
            middle_of_spectrum_in_f1=data.o2p;
        end
        data.of1p=middle_of_spectrum_in_f1;
    end
end
   
if exist([path '/pdata/' num2str(procno) '/proc2'],'file')
    content=textread([path '/pdata/' num2str(procno) '/proc2'],'%s');
    data.phc01_proc2         = str2num(char(content(strmatch('##$PHC0=',content')+1)));
    data.phc11_proc2        = str2num(char(content(strmatch('##$PHC1=',content')+1)));
    data.si1_request        = str2num(char(content(strmatch('##$SI=',content')+1)));
    data.fcor1        = str2num(char(content(strmatch('##$FCOR=',content')+1)));
    
end

if all_quadrants
    run_through_list=0:3;
else
    run_through_list=[0];
end

if exist([path '/pdata/' num2str(procno)  '/2rr'],'file')
    data.typ=1;% not hypercomplex but 2d
end
if exist([path '/pdata/' num2str(procno)  '/2ri'],'file')
    data.typ=2;%  hypercomplex but 2d
end
if exist([path '/pdata/' num2str(procno)  '/2ir'],'file')
    data.typ=2;%  hypercomplex but 2d
end

for run_throug=run_through_list
    switch run_throug
        case 0
            file_to_read='2rr';
        case 1
            file_to_read='2ri';
        case 2
            file_to_read='2ir';
        case 3
            file_to_read='2ii';
       
    end

if exist([path '/pdata/' num2str(procno)  '/' file_to_read],'file')
    if data.ftsize==0
        nbpt1=data.td1;%xf2 data
        % if data.fnmode==tppi_value_fnmode
        %     spectrum=zeros(nbpt1,data.td1/2);
        %  else
        spectrum_in=zeros(nbpt1,data.td1);
        
        % end
    else
        nbpt1=data.si1;%xfb data
        spectrum_in=zeros(nbpt1,si2);
        
    end
    % spectrum=zeros(nbpt1,si2);
    % spectrum=zeros(si2,si1);
    
    file_name_to_read_now=[path '/pdata/' num2str(procno)  '/' file_to_read];
    file_id=fopen(file_name_to_read_now);
    disp(['Reading ' file_name_to_read_now]);
    if data.xdim~=si2
        disp(['Decompressed data (xdim = ' num2str(data.xdim1) ' ' num2str(data.xdim) ])
    end
    %     if data.xdim==si2
    %         for loptd=1:nbpt1
    %             in_r = fread(file_id ,si2,'int');
    %             spectrum(loptd,:)=in_r;
    %         end
    %     else
    
    for l1=0:data.xdim1:(nbpt1-1)
        for l2=1:data.xdim:si2
            for l3=1:2:data.xdim1
                in_1 = fread(file_id ,data.xdim,'int');
                in_2 = fread(file_id ,data.xdim,'int');
                % if (data.ftsize==0) && (data.fnmode==tppi_value_fnmode)% complex data
                %     spectrum((l1+l3-1)/2+1,l2:l2+data.xdim-1)=in_1+1i*in_2;
                %   else
                spectrum_in(l1+l3,l2:l2+data.xdim-1)=in_1;
                spectrum_in(l1+l3+1,l2:l2+data.xdim-1)=in_2;
                % end
            end
        end
    end
    %end
    
    fclose(file_id);
    if data.ftsize==0
         switch run_throug
            case 0
        data.xf2data=spectrum_in/data.nc_proc;
            case 1
                        data.xf2data_i=spectrum_in/data.nc_proc;
            case 2
                        data.xf2data_i=spectrum_in/data.nc_proc;
            case 3
                        data.xf2data_ii=spectrum_in/data.nc_proc;
         end
    else
         switch run_throug
            case 0
                data.spectrum=spectrum_in/data.nc_proc;
            case 1
                data.spectrum_ri=spectrum_in/data.nc_proc;
            case 2
                data.spectrum_ir=spectrum_in/data.nc_proc;
            case 3
                data.spectrum_ii=spectrum_in/data.nc_proc;
         end 
                
       % data.spectrum=spectrum_in/data.nc_proc;
    end
end
end
if exist([path '/pdata/' num2str(procno)  '/1r'],'file')
    file_id=fopen([path '/pdata/' num2str(procno)  '/1r']);
    in_r = fread(file_id ,si2,'int')/nc;
    fclose(file_id);
    data.spectrum=in_r;
    data.typ=3;
end
if exist([path '/pdata/' num2str(procno)  '/1i'],'file')
    file_id=fopen([path '/pdata/' num2str(procno)  '/1i']);
    in_r = fread(file_id ,si2,'int')/nc;
    fclose(file_id);
    data.spectrum_ii=in_r;
        data.typ=3;

end
if isfield(data,'file_not_found')
    file_not_found=data.file_not_found;
end
end