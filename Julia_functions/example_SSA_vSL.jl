# activates the local project environment and downloads all the dependencies.
#ENV["PYTHON"]=""; import Pkg; Pkg.activate(@__DIR__); Pkg.instantiate()
using GpABC, FileIO, JLD2, OrdinaryDiffEq, Distances, Distributions, LinearAlgebra, 
StochasticDiffEq, Statistics, Infiltrator, Plots
pyplot()

include("/Users/rorroagarcia/Desktop/LicensedStatesProject/ABC/Julia/functions.jl")


#structural parameters 
nS=18.3; nL=44.3-nS;
Omega=1;



stoich_mat=[1 -1 1 0; 0 1 -1 -2];
a_vec=collect(vec(1.2:0.2:10.0)); #alpha parameter

#Simulation conditions
time_series=0:0.1:1e3;
initial_conditions=[round(nS),round(nL)];

N_realiz=1000
variancesS_25=zeros(N_realiz,length(a_vec))
variancesSL_25=zeros(N_realiz,length(a_vec))
variancesS_75=zeros(N_realiz,length(a_vec))
variancesSL_75=zeros(N_realiz,length(a_vec))


for (index,value) in enumerate(a_vec)


    N=25;
    println(index)
    @time begin
    a=a_vec[index];
    k1=1; k2=a; k3=(a/(1-nS/N)-1)*nS/nL; 
    k4=nS*(N-nS)/(2*nL*(nL-1));
    reaction_rates=[k1;k2;k3;k4];


    for r in 1:N_realiz
        trajectories=gillespie_algorithm(stoich_mat,
        initial_conditions,reaction_rates,vSP_propensity, time_series,N);

        variancesS_25[r,index]=var(trajectories[:,1]);
        variancesSL_25[r,index]=var(trajectories[:,1].+trajectories[:,2]);

    end

    N=75;    
    k1=1; k2=a; k3=(a/(1-nS/N)-1)*nS/nL; 
    k4=nS*(N-nS)/(2*nL*(nL-1));
    reaction_rates=[k1;k2;k3;k4];


    for r in 1:N_realiz
        trajectories=gillespie_algorithm(stoich_mat,
        initial_conditions,reaction_rates,vSP_propensity, time_series,N);

        variancesS_75[r,index]=var(trajectories[:,1]);
        variancesSL_75[r,index]=var(trajectories[:,1].+trajectories[:,2]);

    end


    end

end

@save "SI_SSA_vSL_2575.jld2" 