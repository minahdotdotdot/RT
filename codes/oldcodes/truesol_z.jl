include("../codes/RT.jl")
seed=1205
Random.seed!(seed);

ϵ = 0.01;        # Nonlinear scale
ω = [-1, 3, -2]; # "slow" wave numbers
C = floatRT(5)   # Energy conserving constants

IC = onUnitCircle(3) # Initial condition
T=1000;               # Final time
L=1000;               # Number of state/amplitudes to save

name="true_seedZ"*string(seed); 
h=exp10(-5);           # time step size
N = Int(ceil(T/h));    # Number of time-steps to get to T
every = Int(ceil(N/L)) # only save L+1 values total
RT_amp(N, h, every, IC, ω=ω, ϵ=ϵ, C=C, stepper=RK4, name=name);