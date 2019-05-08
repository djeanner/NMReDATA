function [first_guess, typical, what_optimizes,min_v, max_v]=get_param_for_sim(optimize,inter,list_diff_chem_shift,list_diff_J1,list_diff_J2,lw_hz,a_h)
disp('Optimization parameters ...')
what_optimizes.amplitude=optimize;
what_optimizes.lw=optimize;
what_optimizes.chemical_shifts=optimize;
what_optimizes.couplings=optimize;
%sum_optimize=what_optimizes.amplitude+what_optimizes.lw+what_optimizes.chemical_shifts+what_optimizes.couplings;

%initialization
first_guess=[];
typical=[];
min_v=[];
max_v=[];

factor=30;
mia=a_h/factor;
maa=factor*a_h;
disp(['Normalization of spectra ' num2str(mia) '-' num2str(maa) ' Hz'])
if what_optimizes.amplitude
    first_guess=[first_guess a_h ];
    typical=[typical a_h];
     min_v=[min_v mia];
    max_v=[max_v maa];
end

milw=0.5;
malw=1.5;
disp(['Line width search ' num2str(milw) '-' num2str(malw) ' Hz'])
if what_optimizes.lw
    first_guess=[first_guess  lw_hz ];
    typical=[typical 1];
    min_v=[min_v milw];
    max_v=[max_v malw];
end

delta_ppm=0.002;
disp(['Chemical shift search +/- ' num2str(delta_ppm) ' ppm'])
if what_optimizes.chemical_shifts
    for lo=1:size(list_diff_chem_shift,1)
        vv=inter.zeeman.scalar{1,list_diff_chem_shift(lo,1)};
        first_guess=[first_guess vv];
        typical=[typical vv];
        min_v=[min_v vv-delta_ppm];
        max_v=[max_v vv+delta_ppm];
    end
    
end

delta_hz=0.01;
disp(['Scalar coupling search +/- ' num2str(delta_hz) ' Hz'])
if what_optimizes.couplings
    for lo=1:size(list_diff_J1,1)
        a=list_diff_J1(lo,1);
        b=list_diff_J2(lo,1);
        
        if a<b
            vv=inter.coupling.scalar{a,b};
        else
            vv=inter.coupling.scalar{b,a};
        end
        first_guess=[first_guess vv];
        typical=[typical vv];
        min_v=[min_v vv-delta_hz];
        max_v=[max_v vv+delta_hz];
    end
end