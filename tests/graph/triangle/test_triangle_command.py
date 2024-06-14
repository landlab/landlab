import os
import subprocess

from landlab.graph.triangle.mesh import TriangleMesh


def test_triangle_noop():
    triangle = TriangleMesh.validate_triangle()

    assert os.path.exists(triangle)

    result = subprocess.run((triangle,), capture_output=True)

    assert result.returncode == 0
    assert result.stdout.startswith(b"triangle")


def test_triangle_example(tmpdir):
    box_poly = """\
# A box with eight vertices in 2D, no attributes, one boundary marker.
8 2 0 1
 # Outer box has these vertices:
 1   0 0   0
 2   0 3   0
 3   3 0   0
 4   3 3   33     # A special marker for this vertex.
 # Inner square has these vertices:
 5   1 1   0
 6   1 2   0
 7   2 1   0
 8   2 2   0
# Five segments with boundary markers.
5 1
 1   1 2   5      # Left side of outer box.
 # Square hole has these segments:
 2   5 7   0
 3   7 8   0
 4   8 6   10
 5   6 5   0
# One hole in the middle of the inner square.
1
 1   1.5 1.5
"""
    triangle = TriangleMesh.validate_triangle()

    with tmpdir.as_cwd():
        with open("box.poly", "w") as fp:
            print(box_poly, file=fp)

        result = subprocess.run(
            (triangle, "box.poly"),
            capture_output=True,
        )

        assert result.stderr == b""
        assert result.stdout.startswith(b"Opening box.poly.")
        assert sorted(os.listdir(".")) == [
            "box.1.ele",
            "box.1.node",
            "box.1.poly",
            "box.poly",
        ]
        assert result.returncode == 0
