# Let's create a 0D energy balance model of Earth's climate
# The simplest model is:
# Change in heat = absorbed stellar radiation - outgoing thermal radiation
# 				+ greenhouse

##
# 1. Solar radiation
# Constants
σ_b = 5.670367 * 10^(-8)  # W m^-2 K^-4
R_sun = 6.957 * 10^(8)  # m
T_eff = 5777.  # K, top of photosphere
distance = 1.496 * 10^(11)  # m
α = 0.3  # reflectivity

luminosity(R, T_eff) = 4 * pi * R^2 * σ_b * T_eff^4;  # [W/m^-2]

irradiance(L, distance) = L / (4 * pi * distance^2)  # [W]

L_sun = luminosity(R_sun, T_eff)
S = irradiance(L_sun, distance) # irradiance 

# Total absorbed solar radiation across planet:
# 'S' hits a disc of area pi*R^2, Earth is ~sphere with area 4*pi*R^2
# So the irradiance is spread out as A_disc / A_Earth = 1/4
# Given an albedo, the amount absorbed scales as (1 - albedo)
absorbed_solar_radiation(; α=α, S=S) = S * (1-α)/4;  # [W / m^2]

# 2. Outgoing thermal radiation
# Consider body cooling to space - G(T) represents the overall
# combined effects of all the complex processes that go into this
# We do a Taylor series expansion to first order:
# G(T) ~= G(T_0) + G'(T_0)(T - T_0) = G'(T_0)T + G(T_0) - T_0*G'(T_0)
# We set a pre-industrial equilibrium temperature (equating incoming and
# outgoing radiation without greenhouse effects) at
T0 = 14.  # [C]
CO2_PI = 280.  # pre-industrial CO2

# To simplify the expression, let's define
# A == G(T_0) - T_0 * G'(T_0)
# B = -G'(T_0)  -> climate feedback parameter
# Hence G(T) ~= A - BT
outgoing_thermal_radiation(T; A=A, B=B) = A - B*T;

# Here we use
B = -1.3  # [W/m^2/C] (negative feedback is stabilising)

# 'A' is found by equating these two terms (absorbed = emitted)
A = S * (1. - α) / 4 + B*T0  # [W/m^2]

# 3. Human greenhouse effect
# Empirically, we know this as a log function of CO2 concentration
greenhouse_effect(CO2; a=a, CO2_PI=CO2_PI) = a * log(CO2 / CO2_PI);

# 4. Change in heat content (LHS of eqn)
# Change in heat capacity (content) is C dT/dt
C = 51.;  # atmosphere and upper-ocean heat capacity [J/m^2/C]


## Numerical solution
# Let's now set up our ODE:
# C dT/dt = (1-α)S/4 - (A-BT) + a ln (CO2/CO2_PI)
# This naturally describes the time evolution of Earth's globally-averaged
# surface temperature
# Now we discretise with finite-differencing the dT/dt term and arrive at
# T_n+1 = T_n + Δt / C [{RHS}]
function timestep!(model)
  append!(model.T, model.T[end] + model.dt * tendency(model));
  append!(model.t, model.t[end] + model.dt);
end

# Re-formulated RHS for time-stepping
tendency(model) = (1. / model.C) * (
  + absorbed_solar_radiation(α=model.α, S=model.S)
  - outgoing_thermal_radiation(model.T[end], A=model.A, B=model.B)
  + greenhouse_effect(model.CO2(model.t[end]), a=model.a, CO2_PI=model.CO2_PI)
);

## Our data structure for the climate model
# We'll use a mutable struct since we want to modify model params
# throughout evolution
mutable struct Model
  T::Array{Float64, 1}

  t::Array{Float64, 1}
  dt::Float64

  CO2::Function

  C::Float64
  a::Float64
  A::Float64
  B::Float64
  CO2_PI::Float64

  α::Float64
  S::Float64
end;

# Let's also make use of multiple-dispatch to define methods that take
# on default values for many parameters
# Constant parameters are optional kwargs
Model(T::Array{Float64, 1}, t::Array{Float64, 1}, dt::Float64, CO2::Function;
      C=C, a=a, A=A, B=B, CO2_PI=CO2_PI, α=α, S=S) = 
(Model(T, t, dt, CO2, C, a, A, B, CO2_PI, α, S));

# Construct from float inputs for convenience
Model(T0::Float64, t0::Float64, dt::Float64, CO2::Function;
      C=C, a=a, A=A, B=B, CO2_PI=CO2_PI, α=α, S=S) = (
        Model([T0], [t0], dt, CO2;
        C=C, a=a, A=A, B=B, CO2_PI=CO2_PI, α=α, S=S)
      );

## Running simulations of the energy balance model
function run!(model::Model, end_year::Real)
  while(model.t[end] < end_year)
    timestep!(model)
  end
end

run!(model) = run!(model, 200.);  # run for 200 years by default

## Let's see what happens with various CO2 functions
CO2_const(t) = CO2_PI;  # constant
