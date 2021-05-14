Multilinear Dynamical System for Tensor Time Series
======

This package provides matlab implementation of [MLDS](https://papers.nips.cc/paper/2013/file/9996535e07258a7bbfd8b132435c5962-Paper.pdf). MLDS is an extension of traditional linear dynamical systems, also known as Kalman filters. It replaces the states and observations to a sequence of tensors instead of vectors. Therefore it is able to handle more complex time series data (e.g. a video clip or graph time series). 

1. to run the example on synthetic data
```
  make demo
```

2. to run mlds on SST data set
```
  matlab -r demo_sst.m
```

Reference
----
```
  Multilinear Dynamical Systems for Tensor Time Series,
  Mark Rogers, Lei Li, and Stuart J. Russell.
  In the 27th Conference on Neural Information Processing Systems(NeurIPS) , 2013.
```
