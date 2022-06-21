### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 8f19895c-78d6-4163-b69d-8373e6cd23d8
@macroexpand @elapsed peakflops()

# ╔═╡ 7208c21c-19f1-4718-9341-8661fbf460d6
begin
	# Functions operate on variables, macros operate on blocks of Julia code
	# This works because Julia code is represetned as Julia objects, which can be
	# manipulated. Then, we can choose what we finally want to evaluate
	x = 1 + 2  # x = 3, the result of evaluating '1 + 2'
	
	# If we quote a value, we can prevent Julia from evaluating it
	expr = :(1 + 2)
	# or
	expr = quote
		1 + 2
	end
	
	typeof(expr)
end

# ╔═╡ Cell order:
# ╠═8f19895c-78d6-4163-b69d-8373e6cd23d8
# ╠═7208c21c-19f1-4718-9341-8661fbf460d6
