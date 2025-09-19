%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LJA 2024
%%%% lewis.alcott@bristol.ac.uk
%%%% model sensitivity analysis

%%%%%% number of runs
sensruns = 100 ;

%%%%%%% multiple runs
parfor N = 1:sensruns
       
    %%%%%%% run model
    fprintf('Run # %d', N )
    fprintf(' of %d \n', sensruns )
    run(N) = Alcott_et_al_2024_NatGeo_front(N) ;

end

%%%%%% define standard time grid for outputs
tgrid = ( run(1).state.time(1) : 1e6 : run(1).state.time(end) ) ;

%%%%%% sens analysis states mapped to tgrid
for N = 1:sensruns
    field_names = fieldnames(run(N).state) ;
    for numfields = 1:length(field_names)
        eval([' sens.' char( field_names(numfields) ) '(:,N) = interp1( run(N).state.time, run(N).state.' char( field_names(numfields) ) ', tgrid) ;'])
    end
end

%%%%%% plotting
CIplots