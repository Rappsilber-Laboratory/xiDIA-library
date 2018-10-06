# xiDIA-library

Tool to generate a cross-linking spectral library for Spectronaut from Xi results.

### Note
This tool requires adaptation in the code regarding the location of the peak list data in order to function in a different environment than the Rappsilber lab.

### Requirements
Python 3

virtualenv

### Installation

Clone git repository to your machine:

```git clone https://github.com/Rappsilber-Laboratory/xiDIA-library.git```

cd into the repository:

```cd xiDIA-library```

Create python virtualenv:

```virtualenv -p python3 --no-site-packages python_env```

Activate virtualenv:

```source python_env/bin/activate```

Install dependencies:

```pip install -r requirements.txt```


### Required files
XiFDR PSM csv file(s)

all runs that are referenced in the XiFDR PSM files in mzML format

### Usage

Modify the config.py.example file, save it as config.py and then run create_spectronaut_lib.py

```
baseDir = "/path/to/your/directory/"    # Path to your working directory
psm_csv_path = baseDir + "psm_csv/"     # Path to subdirectory with your XiFDR psm csv files      
mzml_path = baseDir + "mzml/"           # Path to subdirectory with your mzML files
output_path = baseDir + "out/"          # Path to subdirectory where the library output file will be saved

# iRT equation parameters
# y = m * x + t
# e.g. y= 0.6452 x + 47.469
iRT_m = 0.56                            # iRT equation slope value
iRT_t = 53.34                           # iRT equation y-intercept value

proteinId = "HSA"                       # protein name string
includeNeutralLossFragments = False     # set to True if you want to include neutral loss fragments
writeClSitesToModifiedSequence = True   # set to True to write the cross-linked sites to the ModifiedSequence column
clName = "bs3"                          # cross-linker name string

label_experiment = False                # set to True if you did an heavy-light isotope labeled experiment
clLightLabelName = "bs3-light"          # will be ignored if label_experiment is False
clHeavyLabelName = "bs3-heavy"          # will be ignored if label_experiment is False
deuteriumCount = 4                      # will be ignored if label_experiment is False

fragmentMzLowerLim = 100                # lower m/z limit for fragments
fragmentMzUpperLim = 1000               # upper m/z limit for fragments

# Total number of fragments can not exceed 20!
nClContainingFragments = 15             # maximum number of cross-linker containing fragments that will be included
nLinearFragments = 5                    # maximum number of linear fragments that will be included

```
