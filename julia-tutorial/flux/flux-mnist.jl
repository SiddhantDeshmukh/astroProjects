## MNIST digit recognition with Flux.jl
using Flux, MLDatasets, ImageCore
using Flux.Data: DataLoader
using Statistics
using Flux: crossentropy, throttle, @epochs
using Base.Iterators: repeated

# Load MNIST data
train_x, train_y = MNIST.traindata()
test_x, test_y = MNIST.testdata()

train_loader = DataLoader((train_x, train_y); batchsize = 128, shuffle = true);

# View random image
MNIST.convert2image(MNIST.traintensor(rand(1:60000)))

## Define model
# Simple feed-forward architecture
m = Chain(Dense(784, 40, relu), Dense(40, 10), softmax);

loss(x, y) = crossentropy(m(x), y);
evalcb = () -> @show(loss(test_x, test_y))  # callback to show loss

## Train
@epochs 5 Flux.train!(
  loss,
  params(m),
  train_loader,
  ADAM,
  cb = throttle(evalcb, 10),
);