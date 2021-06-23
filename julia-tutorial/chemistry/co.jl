##
# Let's do a simple CO reaction network to test Catalyst.jl
using Catalyst, DifferentialEquations, Plots, Latexify;
gr();

arnet = @reaction_network begin
	k1, C + O --> CO
	k2, CO + H --> C + O + H
end k1 k2;

function rf(temperature)
	rate = 3e-17 * (temperature ./ 5000) .^ 0.6
	return rate
end

function rd(temperature)
	rate = 1e-10 * (temperature ./ 5000) .^ 22.8
	return rate
end

function abundance_to_number_density(abundance, hydrogen_number_density)
	number_density = hydrogen_number_density * 10^(abundance - 12.)
	return number_density
end

temperature = 6000

# Abundances
A_O = 8.69
A_C = 8.43
A_CO = 8.5

# Parameters [k1, k2]
p = (rf(temperature), rd(temperature))
# p = (1, 2)

# Initial condition [C, O, CO, H]
hydrogen_number_density = 1e12
u₀ = [abundance_to_number_density(A_C, hydrogen_number_density), 
			abundance_to_number_density(A_O, hydrogen_number_density),
 			abundance_to_number_density(A_CO, hydrogen_number_density), 
 			hydrogen_number_density]
# u₀ = [1., 2., 0.5, 5]

# Time interval to solve on
tspan = (1e-6, 1e-2)

# Create ODEProblem
arnetODE = ODEProblem(arnet, u₀, tspan, p)

# Solve and plot
sol = solve(arnetODE, Tsit5())

println("Initial")
println(u₀)
println("Final")
println(last(sol))

plot(sol, lw=2, scale=:log10, legend=:bottomleft)

##
temperatures = LinRange(3000., 8000., 500)
plot(temperatures, rf(temperatures))
plot!(temperatures, rd(temperatures))

##
# Investigate timescales assuming oxygen dominated
function co_eqm_oxygen(temperatures, A_O)
	concentration = (1 .+ rd(temperatures) ./ (rf(temperatures) * A_O)).^-1
	return concentration
end

# Investigate timescales assuming carbon dominated
function co_eqm_carbon(temperatures, A_C)
	concentration = (1 .+ rd(temperatures) ./ (rf(temperatures) * A_C)).^-1
	return concentration
end

# Timescales
function co_timescales(temperatures, A_O, co_eqm, n_H)
	timescales = co_eqm ./ (rf(temperatures) .* A_O .* n_H)
	return timescales
end

A_O = 8e-4
co_eqms_oxygen = co_eqm_oxygen(temperatures, A_O)
co_timescales_oxygen = co_timescales(temperatures, A_O, co_eqms_oxygen, 1.e12)
plot(temperatures, co_eqms_oxygen,
	ylabel="[CO]", xlabel="Temperature [K]")
plot(temperatures, co_timescales_oxygen,
	ylabel="Timescales", xlabel="Temperature [K]")
plot(co_eqms_oxygen, co_timescales_oxygen,
	ylabel="Timescales", xlabel="[CO]")

##