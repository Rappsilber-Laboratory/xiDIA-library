# Paths
baseDir = "/media/data/work/E180125_DDA_Library/"
psm_csv_path = baseDir + "psm_csv/"
mzml_path = baseDir + "mzml/"
output_path = baseDir + "out/"

# iRT equation parameters
# y = m * x + t
# e.g. y= 0.6452 x + 47.469
iRT_m = 0.56
iRT_t = 53.34

proteinId = "HSA"
includeNeutralLossFragments = False
writeClSitesToModifiedSequence = True
clLabelName = "bs3"

label_experiment = False
clLightLabelName = "bs3-light"    # will be ignored if label_experiment is False
clHeavyLabelName = "bs3-heavy"    # will be ignored if label_experiment is False
deuteriumCount = 4              # will be ignored if label_experiment is False

fragmentMzLowerLim = 100
fragmentMzUpperLim = 1000

# Total number of fragments can not exceed 20!
nClContainingFragments = 15
nLinearFragments = 5
