# planar-quad

This is work in progress. On planar local deformations and subdivision approaches.

## Usage

Currently there are three scripts that can be executed.

`$ lein run -m scripts.demo-doo-sabin`

will perform a Doo-Sabin subdivision on a few objects.

`$ lein run -m scripts.demo-catmull-clark`

will perform a Catmull-Clark subdivision on a few objects.

`$ lein run -m scripts.demo-planar-subd`

will perform/attempt a planar subdivision on a few objects.


`$ lein run -m scripts.demo-dual`

will start an interactive jReality session with a 3x3 mesh
that can be dragged on it's boundary. The position of the boundary
will be used to compute the supporting hyperplane for E.

It is suggested that the corner points are not moved first!
During dragging the four points q1 to q4 will computed and in case
of solutions that are convex, have a smaller polygon length than the
boundary and a smaller area than the quadrangle inscribed in the four
corner points the solution is constructed and displayed.

`$ lein run -m scripts.demo-renderer`

will run the render code to render a cube and it's subdivisions
to TikZ code and finally call XeLaTeX on it. It requires Skim and
XeLaTeX as well as the existence of `~/temp`. This can all be
adjusted in `tikz.clj`.

## Installation

It is strongly recommended that [Leiningen](https://github.com/technomancy/leiningen) is used.

`$ lein deps`

`$ lein native-deps`

should be sufficient to install all required dependencies. (JReality,...)

## License

[...]
