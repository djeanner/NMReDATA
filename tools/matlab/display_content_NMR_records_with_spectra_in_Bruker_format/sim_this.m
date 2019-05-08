function [axis_theo_h, spec_theo, spin_system, a_h, lw_hz, inter, parameters, fid_h]=sim_this(sys,inter, bas, list_diff_chem_shift, list_diff_J1, list_diff_J2, dataset, what_optimizes, a_h, lw_hz, param_fit, factor_rescale_to_common_unit)
if nargin>=12
param_fit=param_fit./factor_rescale_to_common_unit;
end

debug=0;

factor_reduce_nb_pt=4;% not very elegant.. should be set according to AQ... not td
% skiped=ones(1,10000);
% apply optimzed params to the system
% lw_hz=param_fit(1,2);
% inter.zeeman.scalar{1,1}=param_fit(1,3);%chemical shift 1

counter=1;
%counter_all=1;
 disp('--------current params---------------')

% if ~skiped(counter_all)
     if what_optimizes.amplitude
        a_h=param_fit(1,counter);
        counter=counter+1;
%         opp=' (optimized)';
%             disp(['scaling factor ' num2str(a_h) opp])

%     else
%         opp='';
%             disp(['scaling factor '  opp])
% 
%     end
    
end
% counter_all=counter_all+1;

% if ~skiped(counter_all)
    if what_optimizes.lw
        lw_hz=param_fit(1,counter);
        counter=counter+1;
%         opp=' (optimized)';
%             disp(['LB ' num2str(lw_hz) ' Hz' opp])
% 
%     else
%         opp='';
%     
%         disp(['LB ' opp])
% 
%     end
    
end
%counter_all=counter_all+1;


 if what_optimizes.chemical_shifts
    for lo=1:size(list_diff_chem_shift,1)
       % if ~skiped(counter_all)
            vv=param_fit(1,counter);
            counter=counter+1;
            for lo2=1:size(list_diff_chem_shift,2)
                indd=list_diff_chem_shift(lo,lo2);
                if indd~=0
                    inter.zeeman.scalar{1,indd}=vv;
                end
            end
%             opp=' (optimized)';
%                     disp(['Chemical shift ' num2str(lo) ' ' num2str(vv) ' ppm' opp])
% 
%         else
%             opp='';
%                     disp(['Chemical shift ' num2str(lo) ' ' opp])
% 
%         end
%                 counter_all=counter_all+1;

    end
end

if what_optimizes.couplings
    for lo=1:size(list_diff_J1,1)
%         if ~skiped(counter_all)
            vv=param_fit(1,counter);
            counter=counter+1;
            for lo2=1:size(list_diff_J1,2)
                for lo3=1:size(list_diff_J2,2)
                    
                    a=list_diff_J1(lo,lo2);
                    b=list_diff_J2(lo,lo3);
                    if (a~=0) && (b~=0)
                        if a<b
                            inter.coupling.scalar{a,b}=vv;
                        else
                            inter.coupling.scalar{b,a}=vv;
                        end
                    end
                end
            end
%             opp=' (optimized)';
%                 disp(['Scalar coupling ' num2str(lo) ' ' num2str(vv) ' Hz' opp])
% 
%         else
%             opp='';
%                 disp(['Scalar coupling ' num2str(lo) ' ' opp])
% 
%         end
% 
%         counter_all=counter_all+1;
%         
     end
end

disp('--------current params---------------')
if what_optimizes.amplitude, opp=' (optimized)'; else opp='';end
disp(['scaling factor ' num2str(a_h) opp])
if what_optimizes.lw, opp=' (optimized)'; else opp='';end
disp(['LB ' num2str(lw_hz) ' Hz' opp])
if what_optimizes.chemical_shifts, opp=' (optimized)'; else opp='';end

for lo=1:size(list_diff_chem_shift,1)
    indd=list_diff_chem_shift(lo,1);
    
    disp(['Chemical shift ' num2str(lo) ' ' num2str(inter.zeeman.scalar{1,indd}) ' ppm' opp])
end
if what_optimizes.couplings, opp=' (optimized)'; else opp='';end
for lo=1:size(list_diff_J1,1)
    a=list_diff_J1(lo,1);
    b=list_diff_J2(lo,1);
    if a<b
        vv=inter.coupling.scalar{a,b};
    else
        vv=inter.coupling.scalar{b,a};
    end
    disp(['Scalar coupling ' num2str(lo) ' ' num2str(vv) ' Hz' opp])
end
if debug

save('./from_sim_this_save_all_variables_before_create.mat')
end
 %   sys.enable={'greedy'};
 % Spinach housekeeping
 spin_system=create(sys,inter);
 if debug

 save('./from_sim_this_save_all_variables_after_create.mat')
 end
 if isfield(inter.coupling,'strength')
     spin_system.inter.coupling.strength=inter.coupling.strength;
 end
 stored_spins_system=spin_system;
 spin_system=basis(spin_system,bas);% error from this when too large spin system
 iso=spin_system.comp.isotopes{1,1};%ake spin of the first spin.....
 warning('change this first may not be the one detected......')
 
 % Sequence parameters - 1H
 parameters.spins={iso};
%               LzSp=state(spin_system,{'Lz','L+'},{1,2});
parameters.rho0=state(spin_system,'L+',iso,'cheap');
parameters.coil=state(spin_system,'L+',iso,'cheap');
%parameters.coil=state(spin_system,'L+','1H','cheap');

parameters.decouple={};
parameters.offseto=dataset.o1p*dataset.sfo1;
% tmp=(dataset.offset-(dataset.sw2/2));
parameters.offset=(dataset.offset-(dataset.sw2/2))*dataset.sfo1;
parameters.sweep=dataset.sw2*dataset.sfo1;
parameters.npoints=dataset.td2/factor_reduce_nb_pt;
parameters.zerofill=dataset.si2;
parameters.axis_units='Hz';
if debug
 save('./from_sim_this_save_all_variables_before_liquid.mat')
end
% Simulation
fid_h=liquid(spin_system,@acquire,parameters,'nmr');

%save('./fid_h.mat','fid_h')

%% Apodization
%dw=1/(2*dataset.sw2h);
%aq=dataset.td2*dw;
%lw_pt=pi*lw_hz*log(dataset.td2/aq)/factor_reduce_nb_pt ;%convert coef  s to pt
lw_pt=pi*lw_hz*log(2*dataset.sw2h)/factor_reduce_nb_pt ;%convert coef  s to pt
% apply line broadening
fid_h=apodization(fid_h,'exp-1d',lw_pt);

% Fourier transform
spec_theo=a_h*real(fftshift(fft(fid_h,parameters.zerofill)));
inc=parameters.sweep/parameters.zerofill;
axis_theo_h=parameters.offset+[-(parameters.sweep-inc)/2 :inc:(parameters.sweep-inc)/2]';%in Hz
axis_theo_h=axis_theo_h/dataset.sfo1;
spec_theo=flipud(spec_theo);
end