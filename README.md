# StrokeStrip
<img src="https://github.com/davepagurek/StrokeStrip/blob/main/img/overview.png?raw=true" />

StrokeStrip jointly parameterizes clusters of strokes (a) that, together, represent strips following a single intended curve (b). We compute the parameterization of this strip (c) restricted to the domain of the input strokes (d), which we then use to produce the parameterized intended curve (d).

## Usage
```sh
./strokestrip input.scap [...args]
```

Additional optional arguments:
 - `--cut`: If your input strokes include sharp back-and-forth turns, this flag will use the Cornucopia library to detect and cut such strokes.
 - `--debug`: Generate extra SVG outputs to introspect the algorithm

## Input format
Drawings are inputted as `.scap` files, which encode strokes as polylines. Strokes are contained in pairs of braces `{ ... }`. Each stroke has a unique stroke id and a cluster id shared by all strokes that colleectively make up one intended curve. Polyline samples can omit pressure by setting it to a default value of 0.

```
#[width]	[height]
@[thickness]
{
	#[stroke_id]	[cluster_id]
	[x1]	[y1]	[pressure1]
	[x2]	[y2]	[pressure2]
	[x3]	[y3]	[pressure3]
	[...etc]
}
[...etc]
```

Example `.scap` inputs are found in the `examples/` directory.

Stroke clusters for new `.scap` files can be generated using the [StrokeAggregator ground truth labeling program.](https://github.com/davepagurek/StrokeAggregatorLabeller)

## Development

### Dependencies
#### Gurobi
This package relies on the Gurobi optimization library, which must be installed and licensed on your machine. If you are at a university, a [free academic license can be obtained.](https://www.gurobi.com/downloads/end-user-license-agreement-academic/) This project was build with Gurobi 9.0; if you are using a newer version of Gurobi, update `FindGUROBI.cmake` to reference your installed version (e.g. change `gurobi90` to `gurobi91` for version 9.1.)

#### Eigen 3
Ensure that Eigen is installed and that its directory is included in `$CMAKE_PREFIX_PATH`.

### Building
StrokeStrip is configured with Cmake:

```sh
mkdir build
cd build
cmake ..
make
```
