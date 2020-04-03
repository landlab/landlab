#! /usr/env/python
"""
Runs RasterModelGrid's unit_test() method
"""

import landlab as ll


def main():
    mg = ll.RasterModelGrid()
    mg._unit_test()


if __name__ == "__main__":
    main()
