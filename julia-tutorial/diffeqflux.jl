##
using DifferentialEquations, DiffEqFlux, Plots

# Lotka-Volterra Equations
# x' = ax - bxy
# y' = -cy + dxy
# Replace nonlinear term 'bxy' with NN
function lotka_volterra!(du, u, p, t)
	x, y = u
	α, β, δ, γ = p
	du[1] = dx = α*x - β * x * y
	du[2] = dx = -δ*x + γ * x * y
end

# Initial condition
u0 = [1.0, 1.0]

# Simulation interval and intermediary points
tspan = (0.0, 10.0)
tsteps = 0.0:0.1:10.0

# Equation parameters
p = [2.0, 1.0, 1.5, 1.25]

# Define problem and solve
prob = ODEProblem(lotka_volterra!, u0, tspan, p)
sol = solve(prob, Tsit5())

# Plot solution
plot(sol)
savefig("LV_ode.png")

# Define loss function for NN
function loss(p)
	sol = solve(prob, Tsit5(), p=p, saveat=tsteps)
	loss = sum(abs2, sol.-1)
	return loss, sol
end

callback = function(p, l, pred)
	display(l)
	plt = plot(pred, ylim=(0, 6))
	display(plt)
	# If 'return true', the optimisation will stop
	# So we return false to tell sciml_train not to halt
	return false
end

result_ode = DiffEqFlux.sciml_train(loss, p, cb=callback, maxiters=100)