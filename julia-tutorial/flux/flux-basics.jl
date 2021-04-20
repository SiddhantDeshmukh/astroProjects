## Basics of Flux.jl
using Plots;
gr();
using Flux;
f(x) = 3x^2 + 2x + 1;
df(x) = gradient(f, x)[1];  # df/dx = 6x + 2
print(df(2));

# Functions with many params
f(x, y) = sum((x .- y) .^ 2);

# Calculate gradients of each param at the same time
gradient(f, [2, 1], [2, 0]);
# For ML models with hundreds of params, work with collections of
# parameters via 'params'
# parameters via 'params'
x = [2, 1];
y = [2, 0];
gs = gradient(params(x, y)) do
  return f(x, y)
end
gs[x];
gs[y];

## Simple linear regression
W = rand(2, 5);
b = rand(2);

predict(x) = W * x .+ b;

function loss(x, y)
  ŷ = predict(x)
  return sum((y .- ŷ) .^ 2)
end

x, y = rand(5), rand(2);  # dummy data
print(loss(x, y));

# Take gradients of W and b wrt loss to improve prediction (grad descent)
gs = gradient(() -> loss(x, y), params(W, b));

# Update 'W' to train model
W̄ = gs[W];
W .-= 0.1 .* W̄
print(loss(x, y));

## Basics of Building Layers
# Let's have two linear layers with a nonlinear sigmoid in between
W1 = rand(3, 5);
b1 = rand(3);
layer1(x) = W1 * x .+ b1;

W2 = rand(2, 3);
b2 = rand(2);
layer2(x) = W2 * x .+ b2;

model(x) = layer2(σ.(layer1(x)));
print(model(rand(5)));  # 2-element vector

## Improved Layer Building
# Lots of repetition here, when we have lots of layers this will be long
# Instead, we can create a function that returns linear layers
function linear(in, out)
  W = randn(out, in)
  b = randn(out)
  return x -> W * x .+ b
end

linear1 = linear(5, 3);  # can access linear1.W, etc
linear2 = linear(3, 2);

model(x) = linear2(σ.(linear1(x)));
print(model(rand(5)));

## Structs
# Equivalently, we could create a struct that explicitly represents
# affine layer
struct Affine
  W::Any
  b::Any
end

Affine(in::Integer, out::Integer) = Affine(randn(out, in), randn(out));

# Overload call so object can be used as a function
(m::Affine)(x) = m.W * x .+ m.b;

a = Affine(10, 5);
print(a(rand(10)));

# And this is basically the 'Dense' layer provided in Flux! It has some
# more useful properties tacked on, but this is its essence

## Chaining Layers
# We can easily chain layers together in Flux
model2 = Chain(Dense(10, 5, σ), Dense(5, 2), softmax)
print(model2(rand(10)));

## Function Composition
# Since models are just functions, we can write the same as a
# function composition
m = Dense(5, 2) ∘ Dense(10, 5, σ);
print(m(rand(10)));

## Chaining Layers
# And 'Chain' works with any Julia function!
m = Chain(x -> x^2, x -> x + 1);
print(m(5));