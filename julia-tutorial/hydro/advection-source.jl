##
# Basic 1D advection
using Plots;
using LinearAlgebra;
gr();

mutable struct LinearAdvection1D
  dataIn::Array{Float64,1}
  N::Int

  Initialize::Function
  CFL::Function
  checkCFL::Function
  upwindMatrixAssembly::Function
  upwindSolve::Function

  function LinearAdvection1D()
    self = new()

    self.Initialize = function (dataIn::Array{Float64,1}, N::Int)
      self.dataIn = dataIn
      self.N = N
    end

    self.CFL = function ()
      deltaX = (self.dataIn[3] - self.dataIn[2]) / self.N
      return abs(self.dataIn[1] * self.dataIn[4] / deltaX)
    end

    self.checkCFL = function ()
      return self.CFL() <= 1 ? true : false
    end

    self.upwindMatrixAssembly = function ()
      alpha_min = min(self.CFL(), 0)
      alpha_max = max(self.CFL(), 0)
      a1 = [alpha_max for n = 1:self.N-1]
      a2 = [1 + alpha_min - alpha_max for n = 1:self.N]
      a3 = [-alpha_min for n = 1:self.N-1]
      A = Tridiagonal(a1, a2, a3) + zeros(self.N, self.N)
      A[1, end] = alpha_max
      A[end, 1] = -alpha_min
      return A
    end

    self.upwindSolve = function (u0::Array{Float64,1})
      return self.upwindMatrixAssembly() * u0
    end

  end
end


# Start of the hydro code
N = 100;  # resolution
x0 = 0;  # left boundary
x1 = 10;  # right boundary
Δt = 0.05;  # time step
c = 1;  # wave speed?
t = 5;  # time

LA1D = LinearAdvection1D()
LA1D.self.Initialize([c; x0; x1; Δt; t], N)
LA1D.self.checkCFL()

x = LinRange(x0, x1, N)
u0 = [exp(-(n - 2.0) * (n - 2.0)) for n in x];

plot(x, u0, title = "Initial Value")

if LA1D.self.checkCFL()
  println("CFL number = $(LA1D.self.CFL())")
  for t = 0:floor(Integer, t / Δt)
    u = LA1D.self.upwindSolve(u0)
    u0 = u
  end

else
  println("CFL number greater than 1! CFL = $(LA1D.self.CFL())")
end

plot(x, u0, title = "Solution at t=$t")
