# Parallelization Tests

## Bezier Surfaces

Results over 20 runs with:
 - `npts = 50`
 - `mpts = 50`
 - `udpts = 2500`
 - `wdpts = 2500`


  computation |   **min**  |   **avg**  |   **max**
--------------|------------|------------|------------
   *serial*   | `3.798914` | `4.449554` | `5.063268`
 *parallel 1* | `4.067732` | `4.496241` | `4.822059`
 *parallel 2* | `4.088791` | `4.522290` | `4.944766`

the test made are:
 - *serial*: no parallelization esplicitly made.
 - *parallel 1*: with `@sync` and `@async` and matrices defined inside `@sync` block.
 - *parallel 2*: with `@sync` and `@async` and matrices defined before `@sync` block.
