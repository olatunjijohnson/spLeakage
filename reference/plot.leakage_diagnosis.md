# Plot a leakage diagnosis

Plot a leakage diagnosis

## Usage

``` r
# S3 method for class 'leakage_diagnosis'
plot(x, which = c("ecdf", "map"), ...)
```

## Arguments

- x:

  A `leakage_diagnosis` object.

- which:

  `"ecdf"` (observed vs target NN-distance ECDFs) or `"map"` (per-point
  leakage contribution in space).

- ...:

  Passed to the underlying plot call.

## Value

`x`, invisibly.
