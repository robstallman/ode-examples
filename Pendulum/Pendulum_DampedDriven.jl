# Set up the package environment
cd(@__DIR__)      # go to the directory of this script
using Pkg         # use the package manager
Pkg.activate(".") # activate the environment defined by the toml files in this directory
Pkg.instantiate() # install missing dependencies, make sure environment is ready to use

# Import packages
using DifferentialEquations # ODEProblem(), solve()
using Plots                 # plot()

# Define parameters
l = 1.0                  # length [m]
m = 1.0                  # mass [kg]
g = 9.81                 # gravitational acceleration [m/s^2]
a = 0.1                  # external torque amplitude [Nm]
Ω = 5.0                  # external torque frequency [rad/s]
M = t -> a * sin(Ω * t)  # external torque [Nm]
b = 1.0                  # damping coefficient [kg * m^2 * s^-1]
p = [l, m, g, M, b] # NOTE: passed into ODEProblem; order persists in the ODEFunction

# Define time range to solve
t_start = 0.0  # [s]
t_end = 10.0   # [s]
tspan = (t_start, t_end) # NOTE: passed into ODEProblem

# Define initial conditions
θ_0 = 0.01  # initial angular deflection [rad]
ω_0 = 0.0   # initial angular velocity [rad/s]
u0 = [θ_0, ω_0] # NOTE: passed into ODEProblem; order persists in the ODEFunction and ODESolution

# Define the function that stores our system of equations.
function pendulum!(du, u, p, t)

    # Rename states for clarity
    θ = u[1]  # [m]
    ω = u[2]  # [m/s]

    # Rename parameters for clarity
    l = p[1]  # length [m]
    m = p[2]  # mass [kg]
    g = p[3]  # gravitational acceleration [m/s^2]
    M = p[4]  # external torque [Nm]
    b = p[5]  # damping coefficient [kg * m^2 * s^-1]

    # Define the differential equations
    #  u[1] = θ, so du[1] = dθ/dt
    #  u[2] = ω, so du[2] = dω/dt = -3g/(2l) sin[θ(t)] - 3b/(ml^2)ω(t) + 3/(ml^2)M(t)
    du[1] = ω
    du[2] = -3g / (2l) * sin(θ) - 3b / (m * l^2) * θ + 3 / (m * l^2) * M(t)

end

# Build the ODE problem
prob = ODEProblem(pendulum!, u0, tspan, p)

# Solve the ODE problem by specifying a numerical method
# NOTE: by passing `save_everystep=true`, it is possible to animate the solution with 
# `animate(solution, idxs=(1,2))` (or any other indexes you're interested in...)
solution = solve(prob, RK4(), save_everystep=true)

# Plot the time-domain trajectory of all state variables on the same set of axes
p1 = plot(solution,  # NOTE: plot() automatically interprets/plots the ODESolution Type
    title="Time Domain: All states variables on one set of axes",
    label=["θ [rad]" "ω [rad/s]"]
)

# Plot the time-domain trajectory of the state variables on subplots
p2 = plot(solution,
    layout=(2, 1),
    label=["θ [rad]" "ω [rad/s]"],
    title=["Time Domain: Angular Deflection" "Time Domain: Angular Velocity"],
    xlabel="Time [s]",
    ylabel=["Angular Deflection" "Angular Velocity"],
    color=["blue" "red"]
)

# Plot the phase portrait of the two state variables
p3 = plot(solution,
    idxs=(1, 2),  # 0 = time, 1 = first state variable, 2 = second state variable, ...
    title="Phase Portrait: Angular Velocity vs Deflection",
    label=["(θ [rad], ω [rad/s])"],
    xlabel="Angular Deflection [rad]",
    ylabel="Angular Velocity [rad/s]"
)