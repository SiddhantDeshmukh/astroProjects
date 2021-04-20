##
# Differential Equations tutorial (Chris Rackauckas)
using DifferentialEquations

# First testing plotting
using Plots; gr()
x = 1:10; y = rand(10);

plot(x, y)
##

# Lotka-Volterra Equations
# x' = ax - bxy
# y' = -cy + dxy
# Define problem
params = (1.0, 2.0, 1.5, 1.25)  # a, b, c, d
f = function(du, u, params, t)  # Define  f as in-place update into du
  a, b, c, d = params
  du[1] = a*u[1] - b*u[1]*u[2]
  du[2] = -c*u[2] + d*u[1]*u[2]
end
u0 = [1.0;1.0]; tspan=(0.0, 10.0)
prob = ODEProblem(f, u0, tspan, params);

# Solve the problem
sol = solve(prob);

# Plot solution using plot recip (Plotly backend)
plot(sol, title="All Plots.jl Attributes are Available")

# Plot phase diagram
plot(sol, title="Phase Diagram", vars=(1,2))

# Let's try some stochastic diffeqs
g = function(du, u, p, t)
  # x' = 0.5x
  # y' = 0.1y
  du[1] = 0.5*u[1]
  du[2] = 0.1*u[2]
end

prob = SDEProblem(f, g, u0, tspan, params)
sol = solve(prob, dt=1/2^4)
plot(sol, title="SDE")
plot(sol, title="Phase Diagram - SDE", vars=(1,2))
