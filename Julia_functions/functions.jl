"""
Runs the Gillespie algorithm for the reaction defined by the stoichiometric matrix, 
reaction rates and propensity functions given, and evaluates the output at the times
given by times_series

Inputs:
- stoichiometric matrix: stoich_mat::Matrix{Int64}
- reaction rates: reaction_rates::Vector{Float64}
- propensity functions: propensity_function()::Vector{Float64}
- evaluation times: time_series::Vector{Float64} or StepRangeLen
- initial conditions: initial_conditions::Vector{Int64}

Output:
- stochastic trajectories: state_series::Matrix{Int64}


""" 


function gillespie_algorithm(stoich_mat, initial_conditions,reaction_rates,propensity_functions, time_series,carrying_capacity)
   
   
    # Initialize variables
    number_of_species=size(stoich_mat)[1]
    #num_reactions = size(stoich_mat, 2)
    state_series = zeros(Int64,length(time_series),number_of_species) 
    state_series[1,:]=initial_conditions[:]
    current_state=initial_conditions
    current_time = 0
    index=2

    # Main loop
    while index <= length(time_series)
        # updates propensities
    
        raw_propensities = reaction_rates.*propensity_functions(current_state,carrying_capacity)
        
        normalisation_factor=sum(raw_propensities)
        
        # Generate two random numbers
        r1, r2 = rand(Float64, 2)
        
        #sets time in which a reaction takes place 
        current_time=current_time-log(r1)/normalisation_factor       
        #if new time exceeds evaluation time, assign current molecule 
        #numbers to the evaluation time and continue to the next step
        if current_time > time_series[index]
            state_series[index,:].=current_state;
            index=index+1
        end
            
        # Choose reaction
        #@infiltrate
        cum_sum_propensities = cumsum(raw_propensities)  
        
        if sum(cum_sum_propensities)==0
            
            try
                state_series[index:end,:].=repeat(current_state',outer=length(time_series)-index+1)
            catch
                #@infiltrate
            end
            index=length(time_series)+1
        else 
            reaction_index = findfirst(cum_sum_propensities .>= r2 * normalisation_factor)
            # Update state
            #@infiltrate
            current_state += stoich_mat[:, reaction_index]
        end
        
    end
    
    return state_series
end



"""
Runs the Gillespie algorithm for the reaction defined by the stoichiometric matrix, 
reaction rates and propensity functions given, and evaluates the output until an extinction 
event takes place

Inputs:
- stoichiometric matrix: stoich_mat::Matrix{Int64}
- reaction rates: reaction_rates::Vector{Float64}
- propensity functions: propensity_function()::Vector{Float64}
- evaluation times: time_series::Vector{Float64} or StepRangeLen
- initial conditions: initial_conditions::Vector{Int64}

Output:
- ectinction time: extinction_time::Float64

""" 

function gillespie_extinction_time(stoich_mat, initial_conditions,reaction_rates,propensity_functions,carrying_capacity)
   
   
    # Initialize variables
    #number_of_species=size(stoich_mat)[1]
    #num_reactions = size(stoich_mat, 2)
    current_state=initial_conditions;
    current_time = 0;
    total_molecules=sum(initial_conditions);

    # Main loop
    while total_molecules != 0 # && current_state!=[0;1] #comment the second condition, only valid for SP_modified

        # updates propensities
        raw_propensities = reaction_rates.*propensity_functions(current_state,carrying_capacity)
        
        normalisation_factor=sum(raw_propensities)

        if normalisation_factor==0
            total_molecules=0;
        else 
            # Generate two random numbers
            r1, r2 = rand(Float64, 2)
        
            #sets time in which a reaction takes place 
            current_time=current_time-log(r1)/normalisation_factor       
            
            # Choose reaction
            cum_sum_propensities = cumsum(raw_propensities)  
         
            reaction_index = findfirst(cum_sum_propensities .>= r2 * normalisation_factor)
            current_state += stoich_mat[:, reaction_index]
            total_molecules=sum(current_state);

        end
           
        #this lines are just for the modified SP----------------------------------------------
       # if current_time>1000
       #     current_time=Inf;
       #     total_molecules=0;
       # end
        #------------------------------------------------------------------------------------
       
        
    end
    
    return current_time
end

"""
This one calculates the first passage time for species species_number reaching target_number

""" 

function gillespie_recovery_time(stoich_mat, initial_conditions,species_number,target_number,reaction_rates,propensity_functions,carrying_capacity)
   
   
    # Initialize variables
    #number_of_species=size(stoich_mat)[1]
    #num_reactions = size(stoich_mat, 2)
    current_state=initial_conditions;
    current_time = 0;
    
    # Main loop
    while current_state[species_number] != target_number 

        # updates propensities
        raw_propensities = reaction_rates.*propensity_functions(current_state,carrying_capacity)
                
        normalisation_factor=sum(raw_propensities)
        
        # Generate two random numbers
        r1, r2 = rand(Float64, 2)
        
        #sets time in which a reaction takes place 
        current_time=current_time-log(r1)/normalisation_factor       
        #if new time exceeds evaluation time, assign current molecule 
        #numbers to the evaluation time and continue to the next step
       # if current_time > time_series[index]
       #     state_series[index,:].=current_state;
       #     index=index+1
       # end
            
        # Choose reaction
        cum_sum_propensities = cumsum(raw_propensities)  
         
        reaction_index = findfirst(cum_sum_propensities .>= r2 * normalisation_factor)
        current_state += stoich_mat[:, reaction_index]
        
    end
    
    return current_time
end

"""
This one calculates the first passage time for the total molecule number to reach target_number

""" 

function gillespie_recovery_time_2(stoich_mat, initial_conditions,target_number,reaction_rates,propensity_functions,carrying_capacity)
   
   
    # Initialize variables
    #number_of_species=size(stoich_mat)[1]
    #num_reactions = size(stoich_mat, 2)
    current_state=initial_conditions;
    current_time = 0;
    
    # Main loop
    while sum(current_state) != target_number 

        # updates propensities
        raw_propensities = reaction_rates.*propensity_functions(current_state,carrying_capacity)
                
        normalisation_factor=sum(raw_propensities)
        
        # Generate two random numbers
        r1, r2 = rand(Float64, 2)
        
        #sets time in which a reaction takes place 
        current_time=current_time-log(r1)/normalisation_factor       
        #if new time exceeds evaluation time, assign current molecule 
        #numbers to the evaluation time and continue to the next step
       # if current_time > time_series[index]
       #     state_series[index,:].=current_state;
       #     index=index+1
       # end
            
        # Choose reaction
        cum_sum_propensities = cumsum(raw_propensities)  
         
        reaction_index = findfirst(cum_sum_propensities .>= r2 * normalisation_factor)
        current_state += stoich_mat[:, reaction_index]
        
    end
    
    return current_time
end


#"""
#Distance function for ABC. Euclidean distance in the space of the 
#first three moments, weighted by their standard deviations using a 
#gaussian approximation.
#"""

#function distance_function_2(x,y)
#    mu1_x=mean(x);
#    mu2_x=moment(x,2);
#    mu3_x=moment(x,3);
#    mu1_y=mean(y);
#    mu2_y=moment(y,2);
#    mu3_y=moment(y,3);

#    distance=sqrt(weight1*(mu1_x-mu1_y)^2+
#    weight2*(mu2_x-mu2_y)^2+
#    weight3*(mu3_x-mu3_y)^2)
#    return distance 
#end

function vBD_propensity(state,N)
    propensity=[state*(N-state[1])/N;state[1]]
    return propensity
end


function vSP_propensity(state,N)
    propensity=[state[1]*(N-state[1])/N;state[1];state[2]*(N-state[1])/N;state[2]*(state[2]-1)/N]
    return propensity
end

function SP_propensity(state,Omega)
    propensity=[state[1];state[1];state[2];state[2]*(state[2]-1)/Omega]
    return propensity
end

function SP_propensity_asym(state,Omega)
    propensity=[state[1];state[1];state[1];state[2];state[2]*(state[2]-1)/Omega]
    return propensity
end

function SP_modified_propensity(state,Omega)
    propensity=[state[1];state[1];state[2]*(state[2]-1)/Omega]
    return propensity
end

function vSP_modified_propensity(state,N)
    propensity=[state[1]*(N-state[1])/N;state[1];state[2]*(state[2]-1)/N]
    return propensity
end

function Poisson_BD_propensity(state,Omega)
    propensity=[1;state[1]]
    return propensity
end
