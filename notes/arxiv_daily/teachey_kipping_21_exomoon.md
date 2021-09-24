# Identifying Potential Exomoon Signals with Convolutional Neural Networks

## Summary

Use CNNs to identify exomoon signals in simulated Kepler light curves.
Feed in light curves as Fnorm = (F - min(F) / ~F - min(F)) - 1,
where F is the array of fluxes and ~F is the median value of the
input segment.

## Questions

- They mention that they use simulated data to avoid observation bias
  encountered, but simulate systems which may be "physically possible
  but rare in nature". Are these systems rare enough in the training
  set as well so that the entire set's probability distribution mimics
  what we already know about star-planet systems?

- They produce an equal number of "moon" and "no-moon" systems:
  doesn't this add a bias in the training set (essentially
  a uniform prior)? They mention that binary classification works
  best when the two choices are balanced in the dataset.

- RNNs / LSTMs? Time-series analysis seems to have gone back to CNNs.

- Generally, ensembles work by "polling the crowd" and allowing
  individual members to specialise. Is there such a specialisation
  discernable, e.g. between short- and long-characteristics in the
  time series?
