# Export animation frames as images for video creation

Exports individual PNG frames that can be assembled into a video using
ffmpeg or similar tools.

## Usage

``` r
export_animation_frames(
  sim,
  output_dir = "frames",
  n_frames = 200,
  width = 800,
  height = 600,
  trail_length = 150
)
```

## Arguments

- sim:

  Object of class "vonkarman_sim".

- output_dir:

  Directory for frame images.

- n_frames:

  Number of frames.

- width, height:

  Image dimensions in pixels.

- trail_length:

  Points in trailing history.

## Value

Vector of output filenames.

## Details

After exporting frames, create video with ffmpeg:
`ffmpeg -framerate 30 -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p output.mp4`

## See also

[`create_animated_attractor()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor.md),
[`create_animated_attractor_accumulate()`](https://robustecologies.github.io/tuRbulence/reference/create_animated_attractor_accumulate.md),
[`save_attractor_html()`](https://robustecologies.github.io/tuRbulence/reference/save_attractor_html.md),
[`create_attractor_3d()`](https://robustecologies.github.io/tuRbulence/reference/create_attractor_3d.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Export frames for video creation
sim <- vonkarman_sim(mu = 0.3, n_steps = 50000)
frames <- export_animation_frames(sim, output_dir = "frames", n_frames = 100)
} # }
```
