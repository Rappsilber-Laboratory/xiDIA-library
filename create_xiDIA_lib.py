import pandas as pd
import pymzml
import numpy as np
import json
import urllib.request
import requests
import datetime
import sys
import copy
import re
import os
import codecs
import ssl
from config import *
from pyteomics import mgf
import ntpath


def get_rt_from_scan_id(scan_id, mzml_reader):
    """
    Get the retention time for a scan from the mzml.

    :param scan_id: Scan id
    :param mzml_reader: Mzml reader object
    :return: Retention time
    """
    try:
        rt = mzml_reader[int(scan_id)]['scan start time']
    except KeyError as e:
        print('could not get rt for scan %s' % scan_id)
        print(e)
        rt = 0
    return rt


def cl_species_id(row):
    unsorted_peps = [
        {"seq": row['PepSeq1'], "linkpos": row['LinkPos1']},
        {"seq": row['PepSeq2'], "linkpos": row['LinkPos2']}
    ]

    sorted_peps = sorted(unsorted_peps, key=lambda k: k['seq'])
    charge = row["match charge"]
    pep1 = sorted_peps[0]['seq']
    pep2 = sorted_peps[1]['seq']
    linkpos1 = sorted_peps[0]['linkpos']
    linkpos2 = sorted_peps[1]['linkpos']

    return '{}_{}-{}_{}:{}'.format(pep1, pep2, linkpos1, linkpos2, charge)


def get_annotation_xidb(psm_row):
    """
    Runs the Xi annotator to get the fragment annotation.
    """
    if psm_row.SearchID > 10000:
        url = 'https://xi3.bio.ed.ac.uk/xiAnnotator/annotate/%s/85160-94827-76653-69142/%s/?peptide=%s&peptide=%s&link=%s&link=%s' % (
            int(psm_row.SearchID), int(psm_row.PSMID), psm_row.PepSeq1, psm_row.PepSeq2, int(psm_row.LinkPos1),
            int(psm_row.LinkPos2))
    else:
        url = 'http://129.215.14.63/xiAnnotator/annotate/%s/85160-94827-76653-69142/%s/?peptide=%s&peptide=%s&link=%s&link=%s' % (
            int(psm_row.SearchID), int(psm_row.PSMID), psm_row.PepSeq1, psm_row.PepSeq2, int(psm_row.LinkPos1),
            int(psm_row.LinkPos2))

    ssl._create_default_https_context = ssl._create_unverified_context
    reader = codecs.getreader("utf-8")
    data = json.load(reader(urllib.request.urlopen(url)))
    return data


def create_json_annotation_request(
        peak_list,
        peptide,
        precursor_charge,
        precursor_mz,
        fragment_types=('peptide', 'b', 'y'),
        fragment_tolerance_ppm=10.0,
        cross_linker="BS3",
        custom_settings=False,
        as_dict=False
):
    """
    Create xiANNOTATOR JSON annotation request.

    :param peak_list: list of lists [[mz, int], ...]
    :param peptide: Peptide dict (pep_seq1, pep_seq2, pep_seq1, link_pos1, link_pos2)
    :param precursor_charge: charge state of the precursor
    :param precursor_mz: m/z of the precursor
    :param fragment_types: ['peptide', 'a', 'b', 'c', 'x', 'y', 'z']
    :param fragment_tolerance_ppm: (float) fragment tolerance in ppm to use
    :param cross_linker: either (str) cross-linker used or (float) cross-linker mod mass
    :param custom_settings: (list) custom_settings
    :param as_dict: returns request as dictionary instead of json
    :return: annotation request
    """

    if type(cross_linker) is float:
        cross_linker_mod_mass = cross_linker
    else:
        cross_linker_mod_masses = {
            'BS3': 138.06807961,
            'DSSO': 158.0037648
        }
        try:
            cross_linker_mod_mass = cross_linker_mod_masses[cross_linker]
        except KeyError:
            raise Exception('unknown cross-linker: %s' % cross_linker)

    mod_mass_dict = {
        'bs3loop': 138.06807,
        'bs3nh2': 155.094619105,
        'bs3oh': 156.0786347,
        'cm': 57.021464,
        'ox': 15.994915,
        'dsso': 158.0037648,
        'dssonh2': 175.030313905,
        'dssooh': 176.0143295,
        'sda-hyd': 100.05243,
        'sda-loop': 82.04186484,
    }

    fragment_tolerance_ppm = str(fragment_tolerance_ppm)

    all_mods = []

    # peptides block
    pep1 = [{'aminoAcid': char, 'Modification': ''} for char in peptide['pep_seq1'] if char.isupper()]
    offset = 1
    for match in re.finditer('[^A-Z]+', peptide['pep_seq1']):
        modification = match.group()
        pep1[match.start() - offset]['Modification'] = modification
        offset += len(modification)
        # add to all mods
        if modification not in all_mods:
            all_mods.append(modification)

    pep2 = [{'aminoAcid': char, 'Modification': ''} for char in peptide['pep_seq2']  if char.isupper()]
    offset = 1
    for match in re.finditer('[^A-Z]+', peptide['pep_seq2']):
        modification = match.group()
        pep2[match.start() - offset]['Modification'] = modification
        offset += len(modification)
        # add to all mods
        if modification not in all_mods:
            all_mods.append(modification)

    peptides_block = [{"sequence": pep1}, {"sequence": pep2}]

    # link sites block
    link_sites_block = [
        {'id': 0, 'peptideId': 0, 'linkSite': peptide['link_pos1'] - 1},
        {'id': 0, 'peptideId': 1, 'linkSite': peptide['link_pos2'] - 1}
    ]

    # peak list
    peak_block = [{"mz": float(x[0]), "intensity": float(x[1])} for x in peak_list]

    # annotation block
    annotation_modifications = [{"aminoAcids": ["*"], "id": mod, "mass": mod_mass_dict[mod]} for mod in all_mods]

    ion_types = [{'type': ion.title() + 'Ion'} for ion in fragment_types]

    annotation_block = {
        "fragmentTolerance": {"tolerance": fragment_tolerance_ppm, "unit": "ppm"},
        "modifications": annotation_modifications,
        "ions": ion_types,
        "crosslinker": {"modMass": cross_linker_mod_mass},
        "precursorCharge": precursor_charge,
        "precursorMZ": precursor_mz,
        "custom": custom_settings
    }

    json_dict = {
        "Peptides": peptides_block,
        "LinkSite": link_sites_block,
        "peaks": peak_block,
        "annotation": annotation_block
    }

    if as_dict:
        return json_dict

    return json.dumps(json_dict)


def get_annotation_json(json_request, annotator_url):
    """
    :param json_request: xiAnnotator json request
    :param annotator_url: URL to xiAnnotator webserver
    :return: annotation response JSON
    """
    headers = {'Content-type': 'application/json', 'Accept': 'application/json'}
    r = requests.post(annotator_url, data=json_request, headers=headers, verify=False)
    response_json = r.json()
    return response_json


def replace_mods(sequence, remove_heavy=True):
    """
    encloses modifications with brackets, e.g. ox -> [ox]
    removes heavy labels from modifications (can be turned off with remove_heavy=False
    """
    if remove_heavy:
        return re.sub('([a-z0-9]+)(-d[0-9]+)?([a-z0-9]+)', r'[\1\3]', sequence)
    else:
        return re.sub('([a-z\-0-9]+)', r'[\1]', sequence)


def strip_sequence(sequence):
    return re.sub('([a-z\-0-9_]+)', '', sequence)


def calculate_iRT_value(row, params):
    if hasattr(row, 'irt_group'):
        m, t = params[int(row.irt_group)]
    else:
        m, t = params[0]
    return (row.rt - t) / m


def check_for_isotope_labeling(crosslinker):
    match = re.search('d([0-9])+', crosslinker)
    if match:
        return int(match.group(1))
    else:
        return 0


def invert_modifications(mod_string, lth=True):
    """
    Inverts the modifications for artificially generated transitions

    Parameter: str, modification string
               bool, indicates light to heavy
    Returns: str, same string but swapping the heavy light mods
    """
    light_mods = ["bs3nh2", "bs3oh"]
    heavy_mods = ["bs3-d4nh2", "bs3-d4oh"]

    if lth:
        # replace light with heavy
        for i, j in zip(light_mods, heavy_mods):
            mod_string = mod_string.replace(i, j)
    else:
        # replace light with heavy
        for i, j in zip(light_mods, heavy_mods):
            mod_string = mod_string.replace(j, i)

    return mod_string


def calc_isotope_frag_mz(mz, charge, cl_containing, label, label_mz_difference):
    """
    Adjusts the mass computation for the un-identified labeled sequence
    """
    if cl_containing:
        if label == 0:
            return mz + label_mz_difference / charge
        else:
            return mz - label_mz_difference / charge
    else:
        return mz


def get_all_mods(mod_seq_list):
    """Gets a list of all modifications in the identifications"""
    matches = []
    for i in mod_seq_list:
        matches.extend(re.findall(".*(\[.*\]).*", i))
    unmatches = np.unique(matches)
    return unmatches


def get_lbl_sequence(peps, lbl, label_both=True):
    """
    Creates the 'LabeledSequence' column for Spectronaut.
    """
    seq1 = replace_mods(peps[0]['seq'], remove_heavy=False)
    seq2 = replace_mods(peps[1]['seq'], remove_heavy=False)

    cl_pos1 = int(peps[0]['pepLinkPos']) - 1
    cl_pos2 = int(peps[1]['pepLinkPos']) - 1

    # get the length of modification label and positions
    modseq1_pos = [i.span() for i in re.finditer('([[a-z\-0-9]+])', seq1)]
    modseq2_pos = [i.span() for i in re.finditer('([[a-z\-0-9]+])', seq2)]

    cl_pos1 = add_position_clsite(modseq1_pos, cl_pos1)
    cl_pos2 = add_position_clsite(modseq2_pos, cl_pos2)

    if label_both:
        lbl_sequence = get_cl_string(cl_pos1, seq1, lbl) + "_" + get_cl_string(cl_pos2, seq2, lbl)
    else:
        lbl_sequence = get_cl_string(cl_pos1, seq1, lbl) + "_" + seq2

    #ToDO: inplement checkAA function again?
    return lbl_sequence


def add_position_clsite(modseq, clpos):
    """
    Adjusts the position of the Cl site with modifications
    """
    # determine the position of the cross-link with modifications in the string
    add_pos = 0
    for modi in modseq:
        # print (cl_pos1, modi[0])
        # -1 is needed since we also subtract it from the cross-link site
        if clpos < modi[0] - 1:
            break
        else:
            add_pos = modi[1] - modi[0]
            # print ("add:", add_pos)
        # add offset to cross-linker position
        clpos = clpos + add_pos
    return clpos


def get_cl_string(clpos, pepseq, lbl):
    """
    Inserts the lbl sequence into the peptide sequence
    """
    return pepseq[:clpos + 1] + "[" + lbl + "]" + pepseq[clpos + 1:]


def check_amino_acid(cl_pos1, modseq1, linkable, entry_template):
    """
    checks if the cl amino acid is a valid one
    """
    if cl_pos1 == 0:
        pass
    else:
        if modseq1[cl_pos1:cl_pos1 + 1] in linkable:
            pass
        else:
            sys.exit("Amino acid not in list of valid AA!\nAA:{}\nCL:{}\nSeq:{}\nStripped:{}".format(
                modseq1[cl_pos1:cl_pos1 + 1], cl_pos1, entry_template["FragmentGroupId"],
                entry_template["StrippedSequence"]))


# checks
if nClContainingFragments + nLinearFragments > 20:
    raise Exception("Total number of fragments should not exceed 20!")

# create PSM DataFrame
print("Reading in PSM csv files...")
psm_files = os.listdir(psm_csv_path)

if 'all_psms_with_rt.csv' in psm_files:
    print('Found all_psms_with_rt.csv in psm_csv folder. Using only this as input!')
    psm_files = ['all_psms_with_rt.csv']

df_list = []
for csv_filename in psm_files:
    if csv_filename.endswith(".csv"):
        print("reading %s..." % csv_filename)
        df_list.append(pd.read_csv(psm_csv_path + csv_filename,
                                   thousands=psm_csv_thousands_sep,
                                   quotechar=psm_csv_quote_char,
                                   sep=psm_csv_col_sep))
    else:
        print("skipped file %s" % csv_filename)
psm_df = pd.concat(df_list, ignore_index=True)

# filter out decoys
psm_df = psm_df[psm_df.isTT]


# for offlineXi read in all MGF files
if offlineXi:
    print('reading in MGF files...')
    # MGF files
    mgf_files = os.listdir(mgf_path)
    mgf_readers = {ntpath.basename(f): mgf.read(mgf_path + f) for f in mgf_files}

# read retention times for onlineXi if they are not in the psm file
# (for offlineXi results the retention times will be read out of the MGF files)
else:
    # Read in retention time
    if 'rt' not in psm_df.columns:
        print("rt column not found in psm file.\nReading in in mzML files...")
        # load mzmls
        mzML_reader_map = {}
        for mzML_filename in os.listdir(mzml_path):
            if mzML_filename.endswith(".mzML"):
                print("reading %s..." % mzML_filename)
                mzML_reader_map[mzML_filename[:-5]] = pymzml.run.Reader(mzml_path + mzML_filename,
                                                                        build_index_from_scratch=True)
            else:
                print("skipped file %s" % mzML_filename)

        # annotate psms with retention time
        print("Getting retention time from mzMLs...")
        psm_df['rt'] = psm_df.apply(
            lambda row: get_rt_from_scan_id(row["scan"], mzML_reader_map[row["run"]]), axis=1)

        # writing csv file with retention times
        psm_df.to_csv(psm_csv_path + 'all_psms_with_rt.csv', index=False)

    # transform rts to iRT -> ToDo: what about the nan columns?
    psm_df['iRT'] = psm_df.apply(lambda row: calculate_iRT_value(row, iRT_params), axis=1)

# create the unique crosslink species id column
psm_df['cl_species_id'] = psm_df.apply(lambda row: cl_species_id(row), axis=1)

# filter to best scoring match per cl_species_id. Unique (peptide pair - link - charge) combination
psm_df.sort_values('Score', inplace=True, ascending=False)
best_scores = []
for index, group in psm_df.groupby('cl_species_id'):
    best_scores.append(group.head(1))
best_df = pd.concat(best_scores)
best_df = best_df.dropna(subset=['PSMID', 'PepSeq1', 'PepSeq2', 'LinkPos1', 'LinkPos2'])


# Create library by looping over the PSMs
lib_list = []
for i, (psm_index, psm) in enumerate(best_df.iterrows()):
    print("{}/{} Done.".format(i, best_df.shape[0]))

    # prepare entry template - information that's common to all fragments of a PSM
    if not label_experiment:
        lbl = clName
    else:
        if check_for_isotope_labeling(psm.Crosslinker) == 0:
            lbl = clLightLabelName
        else:
            lbl = clHeavyLabelName

    # sort peptide sequences and corresponding columns alphabetically by peptide sequence
    unsorted_peps = [
        {
            'seq': psm['PepSeq1'],
            'pepLinkPos': psm['LinkPos1'],
            'protLinkPos': psm['ProteinLinkPos1'],
            'protein': psm['Protein1']
        },
        {
            'seq': psm['PepSeq2'],
            'pepLinkPos': psm['LinkPos2'],
            'protLinkPos': psm['ProteinLinkPos2'],
            'protein': psm['Protein2']
        },
    ]
    sorted_peps = sorted(unsorted_peps, key=lambda k: k['seq'])
    # set flag for switched peptide ids
    switched_pep_ids = psm['PepSeq1'] != sorted_peps[0]['seq']

    # create FragmentGroupId - Unique (peptide pair - link - charge) combination
    fragmentGroupId = psm['cl_species_id']

    # create StrippedSequence - 'linearised' & only amino acids
    strippedSequence = strip_sequence(sorted_peps[0]['seq'] + "_" + sorted_peps[1]['seq'])

    # create LabeledSequence
    labeledSequence = get_lbl_sequence(sorted_peps, lbl, label_both=True)

    # create ModifiedSequence
    if writeClSitesToModifiedSequence:
        modifiedSequence = get_lbl_sequence(sorted_peps, lbl, label_both=True)
    else:
        modifiedSequence = replace_mods(sorted_peps[0]['seq'] + "_" + sorted_peps[1]['seq'])

    # create ProteinId, CrosslinkedResidues and LinkId
    proteinId = f"{sorted_peps[0]['protein']}_{sorted_peps[1]['protein']}"
    crosslinkedResidues = f"{sorted_peps[0]['protLinkPos']}_{sorted_peps[1]['protLinkPos']}"
    linkId = proteinId + '-' + crosslinkedResidues

    # get the scan and annotations if it's offlineXi
    if offlineXi:
        try:
            mgf_reader = mgf_readers[psm['PeakListFileName']]
        except KeyError:
            continue

        scan = mgf_reader[psm['ScanId']]
        rt = scan['params']['rtinseconds'] / 60
        if 'irt_group' in psm.keys():
            irt_g = psm['irt_group']
        else:
            irt_g = 0
        irt = (rt - iRT_params[irt_g][1]) / iRT_params[irt_g][0]

        peaklist = zip(scan['m/z array'], scan['intensity array'])
        peptide = {
            'pep_seq1': psm['PepSeq1'],
            'pep_seq2': psm['PepSeq2'],
            'link_pos1': psm['LinkPos1'],
            'link_pos2': psm['LinkPos2'],
        }
        xiAnn_request = create_json_annotation_request(
            peak_list=peaklist,
            peptide=peptide,
            precursor_charge=psm['Charge'],
            precursor_mz=psm['exp m/z'],
            fragment_types=ion_types,
            fragment_tolerance_ppm=ppm_tol,
            cross_linker=psm['CrosslinkerModMass'],
        )
        xiAnn_json = get_annotation_json(xiAnn_request, annotator_url=annotator_url)
        searchId = 'offline'

    else:
        rt = psm.rt
        irt = psm.iRT
        # get annotations
        xiAnn_json = get_annotation_xidb(psm)
        searchId = psm.SearchID

    fragments = xiAnn_json['fragments']
    clusters = xiAnn_json['clusters']
    entries = []

    entry_template = {
        "linkId": linkId,
        "ProteinId": proteinId,
        "StrippedSequence": strippedSequence,
        "iRT": irt,
        "RT": rt,
        "FragmentGroupId": psm['cl_species_id'],
        "PrecursorCharge": int(psm['match charge']),
        "PrecursorMz": psm['match mass']/psm['match charge'] + 1.00794,
        "ModifiedSequence": modifiedSequence,
        "IsotopeLabel": check_for_isotope_labeling(psm.Crosslinker),
        "scanID": psm['scan'],
        "run": psm.run,
        "searchID": searchId,
        "crosslinkedResidues": crosslinkedResidues,
        "LabeledSequence": labeledSequence
    }

    for fragment in fragments:
        try:
            if fragment['name'] == "P+P":
                continue
            if not includeNeutralLossFragments and re.search('Loss', fragment['type']):
                continue
            for clusterId in fragment['clusterIds']:
                firstPeakId = clusters[clusterId]['firstPeakId']
                re_numbers = re.compile('\d+')
                entry = copy.deepcopy(entry_template)
                entry["FragmentCharge"] = clusters[clusterId]['charge']
                entry["FragmentType"] = fragment['name'][0]
                entry["FragmentNumber"] = re_numbers.findall(fragment['name'])[0]
                if switched_pep_ids:
                    if fragment['peptideId'] == 0:
                        entry['FragmentPepId'] = 1
                    if fragment['peptideId'] == 1:
                        entry['FragmentPepId'] = 0
                else:
                    entry['FragmentPepId'] = fragment['peptideId']
                for cluster in fragment["clusterInfo"]:
                    if cluster['Clusterid'] == clusterId:
                        entry["FragmentMz"] = cluster['calcMZ']
                entry["RelativeFragmentIntensity"] = xiAnn_json['peaks'][firstPeakId]['intensity']
                if 'H20' in fragment['name']:
                    entry["FragmentLossType"] = 'H2O'
                elif 'NH3' in fragment['name']:
                    entry["FragmentLossType"] = 'NH3'
                else:
                    entry["FragmentLossType"] = ''
                if re.search('\+P', fragment['name']):
                    entry["CLContainingFragment"] = True
                else:
                    entry["CLContainingFragment"] = False
                if re.search('Loss', fragment['type']):
                    entry["LossyFragment"] = True
                else:
                    entry["LossyFragment"] = False
                entries.append(pd.Series(entry))
        except Exception as e:
            raise Exception(e, fragment)

    entry_df = pd.DataFrame.from_records(entries)
    if len(entry_df) == 0:
        continue

    # apply filter
    entry_df = entry_df[(entry_df.FragmentMz >= fragmentMzLowerLim) &
                        (entry_df.FragmentMz <= fragmentMzUpperLim)]
    entry_df.drop_duplicates(inplace=True)
    if len(entry_df) == 0:
        continue

    # normalize the relative Intensities to 1
    entry_df.sort_values(by="RelativeFragmentIntensity", inplace=True, ascending=False)
    entry_df.RelativeFragmentIntensity = entry_df.apply(
            lambda row: row.RelativeFragmentIntensity / entry_df.RelativeFragmentIntensity.max(), axis=1)

    lib_list.append(entry_df)

lib_df = pd.concat(lib_list)

filtered_fragments_list = []

for index, group in lib_df.groupby(['FragmentGroupId', 'IsotopeLabel']):
    CLContaining = group[group.CLContainingFragment]
    nonCLContaining = group[~group.CLContainingFragment]
    topGroup = pd.concat([CLContaining.head(nClContainingFragments), nonCLContaining.head(nLinearFragments)])

    filtered_fragments_list.append(topGroup)


lib_df = pd.concat(filtered_fragments_list)

if label_experiment:
    missing_peptides = []
    print("adding missing peptides...")

    deuterium_mass = 2.014102
    hydrogen_mass = 1.00794
    mzdiff = deuteriumCount * (deuterium_mass - hydrogen_mass)

    for i, group in lib_df.groupby('FragmentGroupId'):

        # if only heavy or light exists create the missing peptide
        if len(group.IsotopeLabel.unique()) < 2:

            missing_pep = copy.deepcopy(group)
            missing_pep['FragmentMz'] = missing_pep.apply(
                lambda row: calc_isotope_frag_mz(row.FragmentMz, row.FragmentCharge,
                                                 row.CLContainingFragment, row.IsotopeLabel, mzdiff)
                , axis=1)

            pre_charge = missing_pep["PrecursorCharge"].iloc[0]

            if missing_pep['IsotopeLabel'].mean() == deuteriumCount:
                missing_pep['IsotopeLabel'] = 0
                # count the modifications
                nmods = missing_pep['LabeledSequence'].iloc[0].count("bs3")
                missing_pep['PrecursorMz'] = missing_pep['PrecursorMz'].mean() - (mzdiff * nmods) / pre_charge
                missing_pep['LabeledSequence'] = [invert_modifications(
                    i.replace("[bs3-heavy]", "[bs3-light]"), lth=False)
                                                  for i in missing_pep['LabeledSequence']]

                # changes also the modifications strings...
                # missing_pep['ModifiedSequence'] = [invert_modifications(i, lth=False) for i in
                #                                    missing_pep['ModifiedSequence']]
                missing_pep['Artificial'] = "FromHeavy"

            else:
                nmods = missing_pep['LabeledSequence'].iloc[0].count("bs3")
                # generate the heavy version
                missing_pep['IsotopeLabel'] = 4
                missing_pep['PrecursorMz'] = missing_pep['PrecursorMz'].mean() + (mzdiff * nmods) / pre_charge
                missing_pep['LabeledSequence'] = [invert_modifications(i.replace("[bs3-light]", "[bs3-heavy]"), lth=True)
                                                  for i in missing_pep['LabeledSequence']]
                # changes also the modifications strings...
                # missing_pep['ModifiedSequence'] = [invert_modifications(i, lth=True) for i in
                #                                    missing_pep['ModifiedSequence']]
                missing_pep['Artificial'] = "FromLight"

            # append any missing peptides to add them later to the df
            missing_peptides.append(missing_pep)

    lib_df = pd.concat(missing_peptides + [lib_df])


d = datetime.datetime.now()
d = datetime.date.strftime(d, "%m%d%y")
idstr = "_spectronaut_lbls"

# save the library
print("Writing the library...")
lib_df.to_csv(output_path + proteinId + "_" + "_".join(
        [str(x) for x in lib_df.searchID.unique()]) + "nCL{}_nLin{}_lossy{}_lib{}_{}.csv".format(
        nClContainingFragments, nLinearFragments, includeNeutralLossFragments, idstr, d), index=False)


# write out parameters used for creating this library
with open(output_path + "params.txt", "w") as text_file:
    psm_csv_path = baseDir + "psm_csv/"
    for p in iRT_params:
        text_file.write("iRT_params m: %s\t t: %s\n" % p)
    text_file.write("proteinId: %s\n" % proteinId)
    text_file.write("includeNeutralLossFragments: %s\n" % includeNeutralLossFragments)
    text_file.write("writeClSitesToModifiedSequence: %s\n" % writeClSitesToModifiedSequence)
    text_file.write("clName: %s\n" % clName)
    text_file.write("label_experiment: %s\n" % label_experiment)
    text_file.write("clLightLabelName: %s\n" % clLightLabelName)
    text_file.write("clHeavyLabelName: %s\n" % clHeavyLabelName)
    text_file.write("deuteriumCount: %s\n" % deuteriumCount)
    text_file.write("fragmentMzLowerLim: %s\n" % fragmentMzLowerLim)
    text_file.write("fragmentMzUpperLim: %s\n" % fragmentMzUpperLim)
    text_file.write("nClContainingFragments: %s\n" % nClContainingFragments)
    text_file.write("nLinearFragments: %s\n" % nLinearFragments)
    text_file.write("in files mzML:\n%s\n" % "\n".join(mzml_path))
    text_file.write("in files MGF:\n%s\n" % "\n".join(mgf_path))
    text_file.write("in files PSM CSV:\n%s\n" % "\n".join(psm_csv_path))
    text_file.write("psm_csv_thousands_sep: %s\n" % psm_csv_thousands_sep)
    text_file.write("psm_csv_col_sep: %s\n" % psm_csv_col_sep)
    text_file.write("psm_csv_quote_char: %s\n" % psm_csv_quote_char)
    text_file.write("offlineXi: %s\n" % offlineXi)
    text_file.write("offlineXi: %s\n" % ion_types)
    text_file.write("offlineXi: %s\n" % ppm_tol)
