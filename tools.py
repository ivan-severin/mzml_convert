import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import peakutils
import pyopenms as ms


def plot_spectrum(spectrum):
    # plot every peak in spectrum and annotate with it's m/z
    mz, i = spectrum.get_peaks()
    plt.plot(mz, i)
    # plt.text(mz, i, str(mz))

    # for the title add RT and Precursor m/z if available
    title = ''
    if spectrum.getRT() >= 0:
        title += 'RT: ' + str(spectrum.getRT())
    if len(spectrum.getPrecursors()) >= 1:
        title += '   Precursor m/z: ' + str(spectrum.getPrecursors()[0].getMZ())

    plt.title(title)
    plt.ylabel('intensity')
    plt.xlabel('m/z')
    plt.ylim(bottom=0)

    # plt.show()


def make_centroid(exp=ms.MSExperiment) -> ms.MSExperiment:
    exp_centroid = ms.MSExperiment()
    ms.PeakPickerHiRes().pickExperiment(exp, exp_centroid)
    return exp_centroid


def smoothing(exp=ms.MSExperiment) -> ms.MSExperiment:
    exp = ms.MSExperiment(exp)
    gf = ms.GaussFilter()
    param = gf.getParameters()
    param.setValue("gaussian_width", 1.0)  # needs wider width
    gf.setParameters(param)
    gf.filterExperiment(exp)
    return exp


def baseline(exp=ms.MSExperiment) -> ms.MSExperiment:
    """
        Get Baseline
        """
    res = ms.MSExperiment()
    spectra = exp.getSpectra()
    for spectrum in spectra:
        mz, i = spectrum.get_peaks()
        base_line = peakutils.baseline(i)
        base_line_spectrum = ms.MSSpectrum()
        base_line_spectrum.set_peaks([mz, base_line])
        res.addSpectrum(base_line_spectrum)
        # print("len mz: {}, len baseline: {}".format(mz.shape, baseline.shape))
        # print("len spectra: {}".format(len(res.getSpectra())))

    return res


def peaks(exp, thres=0.02, min_dist=3, min_number=0, mz_start=0., mz_end=float("inf")):
    res = ms.MSExperiment()
    for spectrum in exp.getSpectra():
        # Detect Peaks
        # mz, i = spectrum.get_peaks()
        filtered_mz = []
        filtered_int = []

        for mz, i in zip(*spectrum.get_peaks()):
            # print(mz, i)
            if mz_start < mz < mz_end:
                filtered_mz.append(mz)
                filtered_int.append(i)
        # print("X: {}, Y: {}".format(filtered_mz, filtered_int))
        mz = np.array(filtered_mz)
        i = np.array(filtered_int)
        indexes = peakutils.indexes(i, thres=thres, min_dist=min_dist)
        if len(indexes) >= min_number:
            peaks_spectrum = ms.MSSpectrum()
            peaks_spectrum.set_peaks([mz[indexes], i[indexes]])
            peaks_spectrum.setRT(spectrum.getRT())
            res.addSpectrum(peaks_spectrum)
    return res


def set_ms_level(exp: ms.MSExperiment, ms_level: int) -> ms.MSExperiment:
    spectra = []
    for spectrum in exp.getSpectra():
        spectrum.setMSLevel(ms_level)
        spectra.append(spectrum)
    exp.setSpectra(spectra)
    # for spectrum in exp.getSpectra():
    # print("RT: {}, MS Level: {}".format(spectrum.getRT(), spectrum.getMSLevel()))
    return exp


def remove_precursor(exp):
    pass


def set_ms_level_for_file(mzml_file: str, level: int, output_dir: str) -> None:
    print("Processing: {}".format(mzml_file))
    exp = ms.MSExperiment()
    ms.MzMLFile().load(mzml_file, exp)
    # print("Exp MS Level: {}".format(exp.getMSLevel()))
    exp = set_ms_level_for_file(exp, level)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_file = os.path.join(output_dir, os.path.basename(mzml_file))
    exp = set_ms_level_for_file(exp)
    print("Writing: {}".format(output_file))
    ms.MzMLFile().store(output_file, exp)


def normalize(exp: ms.MSExperiment, method="to_one") -> ms.MSExperiment:
    normalizer = ms.Normalizer()
    param = normalizer.getParameters()
    param.setValue("method", method)
    normalizer.setParameters(param)
    res = ms.MSExperiment(exp)
    normalizer.filterPeakMap(res)
    return res


def do_all(mzml_file, output_dir=None):
    exp = ms.MSExperiment()
    ms.MzMLFile().load(mzml_file, exp)
    exp = set_ms_level(exp, 1)
    print("Spectrum size: {}".format(len(exp.getSpectra())))
    exp_centroid = make_centroid(exp)
    print("Centroid, size: {}".format(len(exp_centroid.getSpectra())))
    # exp_normalized = normalize(exp_centroid)
    mz_start = 560.0
    mz_end = 580.0
    print("Filtering m/z from {} to {}".format(mz_start, mz_end))
    exp_peaks = peaks(exp_centroid, thres=0.3, min_number=2, mz_start=mz_start, mz_end=mz_end)
    peaks_count = sum([len(s.get_peaks()) for s in exp_peaks.getSpectra()])
    print("Found peaks: {}".format(peaks_count))
    rt_s = []
    peaks_mz = []
    peaks_i = []
    for pos, spectrum in enumerate(exp_peaks.getSpectra()):
        mz, i = spectrum.get_peaks()
        for mz_, i_ in zip(mz, i):
            # if len(mz) > 1:
            rt_s.append(spectrum.getRT())
            peaks_mz.append(mz_)
            peaks_i.append(i_)
    d = {"RT,s": rt_s, "m/z": peaks_mz, "Intensity": peaks_i}
    df = pd.DataFrame(data=d)
    output_file = os.path.join(output_dir, os.path.basename(mzml_file) + "_peaks_filtered.csv")
    print("Output: {}".format(output_file))
    df.to_csv(output_file)


def detect_peaks(exp: ms.MSExperiment) -> ms.MSExperiment:
    mass_traces = []
    feature_maps = []

    # Sort spectra by RT. If sort_mz=True also sort each peak in a spectrum by m/z
    # exp.sortSpectra(True)
    mtd = ms.MassTraceDetection()
    mtd_params = mtd.getDefaults()
    mtd_params.setValue(
        "mass_error_ppm", 10.0
    )  # set according to your instrument mass error
    mtd_params.setValue(
        "noise_threshold_int", 400000.0
    )  # adjust to noise level in your data
    mtd.setParameters(mtd_params)
    mtd.run(exp, mass_traces, 0)

    mass_traces_split = []
    epd = ms.ElutionPeakDetection()
    epd_params = epd.getDefaults()
    epd_params.setValue("width_filtering", "fixed")
    epd.setParameters(epd_params)
    epd.detectPeaks(mass_traces, mass_traces_split)

    # elution peak detection
    mass_traces_deconvol = []
    epd = ms.ElutionPeakDetection()
    epd_par = epd.getDefaults()
    epd_par.setValue("width_filtering",
                     "fixed")  # The fixed setting filters out mass traces outside the [min_fwhm: 1.0, max_fwhm: 60.0] interval
    epd.setParameters(epd_par)
    epd.detectPeaks(mass_traces, mass_traces_deconvol)

    # feature detection
    feature_map = ms.FeatureMap()  # output features
    chrom_out = []  # output chromatograms
    ffm = ms.FeatureFindingMetabo()
    # ffm_par = ffm.getDefaults()
    # ffm_par.setValue("remove_single_traces", "true")  # remove mass traces without satellite isotopic traces
    # ffm.setParameters(ffm_par)
    ffm.run(mass_traces_deconvol, feature_map, chrom_out)
    feature_map.setUniqueIds()  # Assigns a new, valid unique id per feature

    # hardcode
    feature_map.setPrimaryMSRunPath([
        "data/covert/KB20221007-1-05.mzML".encode()])  # Sets the file path to the primary MS run (usually the mzML file)
    feature_xml_file = ms.FeatureXMLFile()
    feature_xml_file.store("data/res.featureXML".encode(), feature_map)

    feature_maps.append(feature_map)

    # use as reference for alignment, the file with the largest number of features (works well if you have a pooled QC for example)
    ref_index = feature_maps.index(sorted(feature_maps, key=lambda x: x.size())[-1])
    aligner = ms.MapAlignmentAlgorithmPoseClustering()

    trafos = {}

    # parameter optimization
    aligner_par = aligner.getDefaults()
    aligner_par.setValue("max_num_peaks_considered", -1)  # infinite
    aligner_par.setValue("pairfinder:distance_MZ:max_difference", 10.0)  # Never pair features with larger m/z distance
    aligner_par.setValue("pairfinder:distance_MZ:unit", "ppm")
    aligner.setParameters(aligner_par)
    aligner.setReference(feature_maps[ref_index])

    for feature_map in feature_maps[:ref_index] + feature_maps[ref_index + 1:]:
        trafo = ms.TransformationDescription()  # save the transformed data points
        aligner.align(feature_map, trafo)
        trafos[feature_map.getMetaValue("spectra_data")[0].decode()] = trafo
        transformer = ms.MapAlignmentTransformer()
        transformer.transformRetentionTimes(feature_map, trafo, True)

    # # align mzML files based on FeatureMap alignment and store as mzML files (for GNPS!)
    # for file in mzML_files:
    #     exp = ms.MSExperiment()
    #     ms.MzMLFile().load(file, exp)
    #     exp.sortSpectra(True)
    #     exp.setMetaValue("mzML_path", file)
    #     if file not in trafos.keys():
    #         ms.MzMLFile().store(file[:-5]+"_aligned.mzML", exp)
    #         continue
    #     transformer = ms.MapAlignmentTransformer()
    #     trafo_description = trafos[file]
    #     transformer.transformRetentionTimes(exp, trafo_description, True)
    #     ms.MzMLFile().store(file[:-5]+"_aligned.mzML", exp)
    # mzML_files = [file[:-5]+"_aligned.mzML" for file in mzML_files]
