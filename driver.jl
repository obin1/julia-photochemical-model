#################### Part 1: generate reference data ####################
using DiffEqBase, OrdinaryDiffEq
using Catalyst
using Random
using JLD
using Surrogates
export ChemSys, random_sims
include("jpm-mech.jl")

ref_model = JuliaPhotochem.ChemSys()
nspecs = 11

equations(convert(ODESystem, ref_model.rn))

ndays = 3
timelength = 60 * (ndays * 24) # minutes
saveat = 60.0 # minutes
JuliaPhotochem.random_sims(ref_model, 3000, timelength, saveat, 1:3000)


#     @time ref_data_spinup_train, ref_params_train = JuliaPhotochem.random_sims(ref_model, 3000, timelength, saveat, 1:3000);
#     timesteps_keep_train = ref_data_spinup_train[1].t .> 24*60
#     times_train = ref_data_spinup_train[1].t[timesteps_keep_train]
#     time_train = ref_data_spinup_train[1].t
#     want_specs = setdiff(1:13, [ref_model.specmap["O2"], ref_model.specmap["H2O"]])
#     for j in 1:10
#         for i in 1:length(ref_data_spinup_train)
#             if ref_data_spinup_train[i].retcode != :Success
#                 ref_data_spinup_train[i] = ref_data_spinup_train[i-1]
#             end
#         end
#     end

#     ref_data_train = ref_data_spinup_train[want_specs, timesteps_keep_train, :]

    
#     JLD2.jldsave(ProjectPath*"ref_data_train.jld";ref_data_train, ref_params_train, times_train)
#    GC.gc()


#    ndays = 10
#    timelength = 60 * (ndays * 24) # minutes
#    saveat = 60.0 # minutes
#    @time ref_data_spinup_validate, ref_params_validate = random_sims(ref_model, 375, timelength, saveat, 3001:3375);
#    timesteps_keep_validate = ref_data_spinup_validate[1].t .> 24*60
#    times_validate = ref_data_spinup_validate[1].t[timesteps_keep_validate]
#    time_validate = ref_data_spinup_validate[1].t
#    for j in 1:10
#        for i in 1:length(ref_data_spinup_validate)
#            if ref_data_spinup_validate[i].retcode != :Success
#                ref_data_spinup_validate[i] = ref_data_spinup_validate[i-1]
#            end
#        end
#    end

#    ref_data_validate = ref_data_spinup_validate[want_specs, timesteps_keep_validate, :]

   
#    JLD2.jldsave(ProjectPath*"ref_data_validate.jld";ref_data_validate, ref_params_validate, times_validate)
#   GC.gc()





#   @time ref_data_spinup_test, ref_params_test = random_sims(ref_model, 375, timelength, saveat, 3376:3750);
#   timesteps_keep_test= ref_data_spinup_test[1].t .> 24*60
#   times_test = ref_data_spinup_test[1].t[timesteps_keep_test]
#   time_test = ref_data_spinup_test[1].t
#   for j in 1:10
#       for i in 1:length(ref_data_spinup_test)
#           if ref_data_spinup_test[i].retcode != :Success
#               ref_data_spinup_test[i] = ref_data_spinup_test[i-1]
#           end
#       end
#   end

#   ref_data_test = ref_data_spinup_test[want_specs, timesteps_keep_test, :]
