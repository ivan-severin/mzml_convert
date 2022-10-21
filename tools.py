import peakutils
import pyopenms as ms
import matplotlib.pyplot as plt


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


def peaks(exp, thres=0.3, min_dist=5, window=15, order=5, limit_peaks=False):
    res = ms.MSExperiment()
    for spectrum in exp.getSpectra():
        # Detect Peaks
        mz, i = spectrum.get_peaks()
        indexes = peakutils.indexes(i, thres=thres, min_dist=min_dist)
        # print(mz[indexes], i[indexes])
        peaks_spectrum = ms.MSSpectrum()
        peaks_spectrum.set_peaks([mz[indexes], i[indexes]])
        res.addSpectrum(peaks_spectrum)
    return res


def do_all(mzml_file):
    exp = ms.MSExperiment()
    ms.MzMLFile().load(mzml_file, exp)
    # exp_centroid = make_centroid(exp)
    exp_centroid = exp
    # exp_baseline = baseline(exp_centroid)
    exp_peaks = peaks(exp_centroid, thres=0.3)
    return exp_centroid, exp_peaks
    #
    # if limit_peaks and len(self.indexes) > 4:
    #     self.indexes = peakutils.indexes(self.y_smooth, thres=thres + 0.3, min_dist=min_dist + 20)

# def detect_peaks(exp=ms.MSExperiment) -> ms.MSExperiment:
#     mass_traces = []
#
#     # Sort spectra by RT. If sort_mz=True also sort each peak in a spectrum by m/z
#     # exp.sortSpectra(True)
#     mtd = ms.MassTraceDetection()
#     mtd_params = mtd.getDefaults()
#     mtd_params.setValue(
#         "mass_error_ppm", 10.0
#     )  # set according to your instrument mass error
#     mtd_params.setValue(
#         "noise_threshold_int", 4000.0
#     )  # adjust to noise level in your data
#     mtd.setParameters(mtd_params)
#     mtd.run(exp, mass_traces, 0)
#
#     mass_traces_split = []
#     mass_traces_final = []
#     epd = ms.ElutionPeakDetection()
#     epd_params = epd.getDefaults()
#     epd_params.setValue("width_filtering", "fixed")
#     epd.setParameters(epd_params)
#     epd.detectPeaks(mass_traces, mass_traces_split)
#
#     # elution peak detection
#     mass_traces_deconvol = []
#     epd = ms.ElutionPeakDetection()
#     epd_par = epd.getDefaults()
#     epd_par.setValue("width_filtering",
#                      "fixed")  # The fixed setting filters out mass traces outside the [min_fwhm: 1.0, max_fwhm: 60.0] interval
#     epd.setParameters(epd_par)
#     epd.detectPeaks(mass_traces, mass_traces_deconvol)
#
#     # feature detection
#     feature_map = ms.FeatureMap()  # output features
#     chrom_out = []  # output chromatograms
#     ffm = ms.FeatureFindingMetabo()
#     ffm_par = ffm.getDefaults()
#     ffm_par.setValue("remove_single_traces", "true")  # remove mass traces without satellite isotopic traces
#     ffm.setParameters(ffm_par)
#     ffm.run(mass_traces_deconvol, feature_map, chrom_out)
#     feature_map.setUniqueIds()  # Assigns a new, valid unique id per feature
#     # feature_map.setPrimaryMSRunPath([file.encode()]) # Sets the file path to the primary MS run (usually the mzML file)
#     # feature_maps.append(feature_map)
