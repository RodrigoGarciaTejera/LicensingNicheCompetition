%%-------------------------------------------------------------------------
% Rodrigo García-Tejera
% Contact info: rodrigo.garcia@ed.ac.uk ; rgarcia@fisica.edu.uy
% Affiliations: 
% Centre for Regenerative Medicine, University of Edinburgh, UK. 
% Created: 19-12-2022
% Updated: 20-12-2022
%%-------------------------------------------------------------------------


classdef RxnNet 
    
    properties
        stoichiometricMatrix
        rateVector
        propensityFunctions
    end

    methods

    %Constructor for BirthDeathProcess objects
    function obj= RxnNet(stoichiometricMatrix,rateVector,propensityFunctions)
            
            obj.stoichiometricMatrix=stoichiometricMatrix;
            obj.rateVector=rateVector;
            obj.propensityFunctions=propensityFunctions;
        
     end     
        
     function nSpecies=numberOfSpecies(obj)
        nSpecies=numel(obj.stoichiometricMatrix(:,1));    
     end



%Runs the stochastic simulation algorithm (SSA) for the chemical reaction 
%network defined in obj [1]. 
%
%Inputs:  
% obj - chemical reaction network object
% timeVector - vector of times where the molecule numbers are going to be
%              evaluated
% initialNumbers - vector of initial molecule numbers
% extinction - boolean variable to state whether there is an absorbing
% boundary at the zero molecule state. To speed up the simulation if 
% extinction=true, simulation ends whenever the sum of all molecule numbers
% is zero.  
%
% Ouput:
% trajectories - matrix with stochastic trajectories, one column per 
% species. size(trajectories) = numel(timeVector) x number of species.
%
% Note: Beware of the sampling rate introduced by timeVector. Oversampling
% can result in a slight bias.
%
% [1] Gillespie, Daniel T. "Exact stochastic simulation of coupled chemical
% reactions." The journal of physical chemistry 81.25 (1977): 2340-2361.
  

function trajectories=SSA(obj,timeVector,initialNumbers,extinction,varargin)
            

        if nargin~=0
            N=varargin;
        end

        %number of species in the reaction network
        [numberOfSpecies,~]=size(obj.stoichiometricMatrix);

        
        %memory allocation for the stochastic trajectories
        trajectories=nan(numel(timeVector),numberOfSpecies);
        
        %initial conditions
        moleculeNumbers=initialNumbers;
        
        currentTime=0;
        
        %--------------ITERATION OF THE GILLESPIE ALGORITHM----------------
        
        %index for the time vector
        index=1; 
    
        while index<=numel(timeVector)

            %updates propensities  
            try
            rawPropensities=obj.rateVector.*obj.propensityFunctions(moleculeNumbers);
            catch
                keyboard
            end
                             
            %updates normalisation factor
            propensityNormalisation=sum(rawPropensities);                     

            %picks two random numbers uniformly distributed in [0,1] 
            draw=rand(1,2);    
                
            %sets time in which a reaction takes place 
            currentTime=currentTime-log(draw(1))/propensityNormalisation;
            
            %if new time exceeds evaluation time, assign current molecule 
            % numbers to the evaluation time and continue to the next step
            if currentTime > timeVector(index)
               trajectories(index,:)=moleculeNumbers;
               index=index+1;
            end
            
            %chooses which reaction is going to take place
            cumProbabilities=cumsum(rawPropensities/propensityNormalisation);
            I=find(cumProbabilities-draw(2)>0,1,'first'); 

            %if isempty(I)
            %    keyboard
            %end
            
            %updates molecule numbers
            moleculeNumbers=moleculeNumbers+obj.stoichiometricMatrix(:,I)';
            
            
            %the following is just to speed up the simulation when n=0 is an
            %absorbing boundary. 
            if extinction==true
                if sum(moleculeNumbers)==0 
                    trajectories(index:end,:)=0;
                    index=numel(timeVector)+1;
                end
            end
        end
    end

%Calculates numerically the size of fluctuations according to the linear 
%noise approxmation for the reaction network object and the Jacobian in 
%steady state, assuming that a steady state exists. Calculation of LNA is
%done solving Lyapunov equation, see Elf, J., & Ehrenberg, M. (2003). Fast 
% evaluation of fluctuations in biochemical networks with the linear noise 
% approximation. Genome research, 13(11), 2475-2484.

%Inputs:  - obj: RxNet object
%         - steadyStates: (vector) steady state in which the LNA is assesed
%         - Omega: System volume
%Output:  - fluctuationsMatrix: covariance matrix  


function fluctuationsMatrix=LNAnumerical(obj,steadyStates,Omega,propensityFunctionsDet)

    %Construction of the Jacobian matrix in steady state-------------------
    
    %creating symbolic variables for each species
    nSpecies=obj.numberOfSpecies;%number of species in the reaction network
    x=[];
    for k=1:nSpecies
        eval(['syms phi',num2str(k),';']);
        eval(['x=[x,phi',num2str(k),'];']);
    end
    
    %deterministic propensity functions
    f=propensityFunctionsDet(x);
   
    %Build up of the rate equations
    F=[];
    for k=1:nSpecies
        F=[F;f.*obj.stoichiometricMatrix(k,:).*obj.rateVector];
    end
    rateEqs=sum(F,2);%symbolic expression of the rate equations' transfer 
    %matrix
    
    %calculation of the Jacobian
    Jacobian=[];%Jacobian of the rate equations
    for k=1:nSpecies
        eval(['Jacobian=[Jacobian,diff(rateEqs,phi',num2str(k),')];']);%differentiates 
        % against every variable
    end
    for k=1:nSpecies
        eval(['phi',num2str(k),'=steadyStates(k);'])
    end
    Jacobian=eval(Jacobian); %this is the final expression for the Jacobian 
    %in steady state

    
    %I COULD CALCULATE THE STEADY STATES SYMBOLICALLY, BUT IT WOULDNT
    %WORK IF THE RATE EQUATIONS ARE MULTISTABLE, THAT'S WHY I ENTER THE 
    %STEADY STATES AS PARAMETERS.

    %NOTE THAT THE STOCHASTIC AND DETERMINIST PROPENSITY FUNCTIONS DIFFER
    %WHEN THERE IS A MULTIMOLECULAR REACTION THAT INVOLVES MORE THAN ONE
    %MOLECULE OF THE SAME SPECIES, LIKE 2P -->0 IF THE STEADY STATE VALUES 
    %IN THAT VARIABLE IS BIG ENOUGH IT SHOULDNT CHANGE MUCH. 
       
    %----------------------------------------------------------------------

    % Diffusion matrix
    diffusionMatrix=Omega*obj.stoichiometricMatrix*diag(eval(f.*obj.rateVector)) ... 
    *obj.stoichiometricMatrix';%diffusion matrix
    % Solution to Lyapunov equation
    fluctuationsMatrix=lyap(Jacobian,diffusionMatrix); %uses lyapunov solver 
    %from the control systems toolbox
    if fluctuationsMatrix(1,1)<0
        keyboard
    end

    end              
    

    end
end
