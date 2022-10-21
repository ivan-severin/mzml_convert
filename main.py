#!/usr/bin/env python

import sys

import numpy as np

from tools import *
import pandas as pd

def get_integral_spectrum(exp) :
    # res = np.ndarray()
    for i in range(exp.getSize()):
        print(i)



def main(mzml_file, output_file):
    """
    usage:

        ./simple_parser.py <path_to_mzml_file> <output_csv_file>

    """
    exp_centroid, exp_peaks = do_all(mzml_file)

    rt_s = [spectrum.getRT() for spectrum in exp_centroid.getSpectra()]
    peaks_mz = [spectrum.get_peaks()[0] for spectrum in exp_peaks.getSpectra()]
    indexes = np.arange(0, exp_centroid.size())
    print("RTs size: {}".format(len(rt_s)))
    print("Peaks size: {}".format(len(peaks_mz)))
    print("Indexes size: {}".format(len(indexes)))
    d = {'RT,s': rt_s, "indexes": indexes, "peaks_mz": peaks_mz}
    df = pd.DataFrame(data=d)
    df.to_csv(output_file)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(main.__doc__)
        exit()
    mzml_file = sys.argv[1]
    output_file = sys.argv[2]
    main(mzml_file, output_file)
