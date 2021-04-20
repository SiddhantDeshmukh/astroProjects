## Recurrent Cells
using Flux
# From a simple feed-forward y_1 = f(x_1), let's now implement a recurrent
# y_1 = f(h, x_1) where 'h' is a hidden state
# Flux's support closely follows the mathematical perspective
Wxh = randn(5, 10);
Whh = randn(5, 5);
b = randn(5);

function rnn(h, x)
  h = tanh.(Wxh * x .+ Whh * h .+ b)
  return h, h
end

x = rand(10);  # dummy data
h = rand(5);  # initial hidden state
h, y = rnn(h, x);
print(h, y);

## Flux RNN
# Equivalently, the RNN in flux is called with
rnn2 = Flux.RNNCell(10, 5);

x = rand(10);  # dummy data
h = rand(5);  # initial hidden state

h, y = rnn2(h, x);
print(h, y);

## Stateful Models
# Use the 'Recur' wrapper to treat hidden layers as stateful
m = Flux.Recur(rnn, h);
y = m(x);
print(y);

## RNN constructor
# This is just a wrapped cell!
print(RNN(10, 5));

## Sequences
# How do we work with sequences of data? We can use 'Recur' to apply our
# model to each element of a sequence!
seq = [rand(10) for i = 1:10];
m.(seq)  # returns a list of 5-element vectors