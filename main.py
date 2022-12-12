#!/usr/bin/env python

import sys

import numpy as np

from tools import *
import pandas as pd

def get_integral_spectrum(exp):
    # res = np.ndarray()
    for i in range(exp.getSize()):
        print(i)


def main(mzml_file, output_dir):
    """
    usage:

        ./main.py <mzml_files...> <output_dir>

    """
    # exp_centroid, exp_peaks = do_all(mzml_file)

    do_all(mzml_file, output_dir)




if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
        exit()
    mzml_file = sys.argv[1]
    output_dir = sys.argv[2]
    main(mzml_file, output_dir)
