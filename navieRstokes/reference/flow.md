# Create animation from navieRstokes simulation results

S3 method to generate an MP4 video or GIF animation from Navier-Stokes
simulation results. Frames are rendered through the ggplot2-based
`plot.navieRstokes` method and assembled by either the av (MP4) or
gifski (GIF) backend.

## Usage

``` r
flow(x, ...)

# S3 method for class 'navieRstokes'
flow(
  x,
  time_seq = NULL,
  output_file = "flow_animation.mp4",
  fps = 10,
  format = "mp4",
  plot_type = "vorticity",
  width = 800,
  height = 600,
  res = 100,
  ...
)
```

## Arguments

- x:

  Object of class 'navieRstokes' (output from
  [`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)).

- ...:

  Passed to method.

- time_seq:

  Integer vector of time indices to include. Defaults to every saved
  snapshot.

- output_file:

  Character. Output filename (e.g. `"flow.mp4"`, `"flow.gif"`).

- fps:

  Numeric. Frames per second.

- format:

  `"mp4"` or `"gif"`.

- plot_type:

  Forwarded to `plot.navieRstokes`; one of `"vorticity"`, `"pressure"`,
  `"speed"`, `"velocity"`.

- width, height, res:

  Frame width, height and resolution in pixels.

## Value

Invisibly returns the path to the output file.

## See also

`flow.navieRstokes`,
[`plot.navieRstokes`](https://robustecologies.github.io/navieRstokes/reference/plot.navieRstokes.md),
[`simulate_navier_stokes`](https://robustecologies.github.io/navieRstokes/reference/simulate_navier_stokes.md)
