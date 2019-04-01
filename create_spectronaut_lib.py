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


def get_rt_from_scan_id(scan_id, mzml_reader):
    try:
        # rt = msrun[id1]['scan start time']
        rt = mzml_reader[int(scan_id)]['scan start time']
    except KeyError as e:
        print('could not get rt for scan %s' % scan_id)
        print(e)
        rt = 0
    return rt


def create_unique_pep_seq(row):
    unsorted_peps = [
        {"seq": row['PepSeq1'], "linkpos": row['LinkPos1']},
        {"seq": row['PepSeq2'], "linkpos": row['LinkPos2']}
    ]

    sorted_peps = sorted(unsorted_peps, key=lambda k: k['seq'])
    charge = row["match charge"]
    # crossLinker = row['Crosslinker']
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
    Create XiAnnotator JSON annotation request
    :param peak_list: list of lists [[mz, int], ...]
    :param peptide: Peptide Object
    :param precursor_charge:
    :param precursor_mz:
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
            'DSSO': 158.0037648,
            "UCL2": 250.012224,
            "UCL3": 218.040152
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
        'dssooh': 176.0143295
    }

    fragment_tolerance_ppm = str(fragment_tolerance_ppm)

    all_mods = []

    # peptides block
    pep1 = [{'aminoAcid': char, 'Modification': ''} for char in peptide.pep_seq1 if char.isupper()]
    offset = 1
    for match in re.finditer('[^A-Z]+', peptide.pep_seq1):
        modification = match.group()
        pep1[match.start() - offset]['Modification'] = modification
        offset += len(modification)
        # add to all mods
        if modification not in all_mods:
            all_mods.append(modification)

    pep2 = [{'aminoAcid': char, 'Modification': ''} for char in peptide.pep_seq2 if char.isupper()]
    offset = 1
    for match in re.finditer('[^A-Z]+', peptide.pep_seq2):
        modification = match.group()
        pep2[match.start() - offset]['Modification'] = modification
        offset += len(modification)
        # add to all mods
        if modification not in all_mods:
            all_mods.append(modification)

    peptides_block = [{"sequence": pep1}, {"sequence": pep2}]

    # link sites block
    link_sites_block = [
        {'id': 0, 'peptideId': 0, 'linkSite': peptide.link_pos1 - 1},
        {'id': 0, 'peptideId': 1, 'linkSite': peptide.link_pos2 - 1}
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


def get_annotation_json(json_request):
    """
    :param json_request: xiAnnotator json request
    :return: annotation response JSON
    """
    headers = {'Content-type': 'application/json', 'Accept': 'application/json'}
    url = 'https://xi3.bio.ed.ac.uk/xiAnnotator/annotate/FULL'
    r = requests.post(url, data=json_request, headers=headers)
    response_json = r.json()
    return response_json


# def replace_HL_mods(sequence):
#	 """ Removes heavy/light modifications of the cross-linker, e.g. bs3nh2 """
#	 mods = ["bs3nh2", "bs3oh", "bs3-d4nh2", "bs3-d4oh"]
#	 mods = "|".join(mods)
#	 return replace_mods(re.sub('({})'.format(mods), '', sequence))


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


def calculate_iRT_value(rt, m, t):
    return (rt - t) / m


def check_for_isotope_labeling(crosslinker):
    match = re.search('d([0-9])+', crosslinker)
    if match:
        return int(match.group(1))
    else:
        return 0


# def sanity_check_RT():
#	 expmz = []
#	 scanmz = []
#	 for ii, row in psm_df.iterrows():
#
#		 try:
#			 scanmz.append(float(msruns[row["run"]][row["scan"]]["precursors"][0]["mz"]))
#			 expmz.append(row["exp m/z"])
#		 except:
#			 pass
#
#	 expmz = np.array(expmz)
#	 scanmz = np.array(scanmz)
#	 diff = scanmz - expmz
#	 plt.hist(diff)


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


def get_lbl_sequence(psm, lbl, label_both=True):
    """
    Creates the 'LabeledSequence' column for spectronaut.

    """

    seq1 = replace_mods(psm.PepSeq1, remove_heavy=False)
    seq2 = replace_mods(psm.PepSeq2, remove_heavy=False)

    cl_pos1 = int(psm.LinkPos1) - 1
    cl_pos2 = int(psm.LinkPos2) - 1

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
df_list = []
for csv_filename in os.listdir(psm_csv_path):
    if csv_filename.endswith(".csv"):
        print("reading %s..." % csv_filename)
        df_list.append(pd.read_csv(psm_csv_path + csv_filename))
    else:
        print("skipped file %s" % csv_filename)
psm_df = pd.concat(df_list, ignore_index=True)

# load mzmls
print("Reading in mzML files...")
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

# transform rts to iRT -> what about the nan columns?
psm_df['iRT'] = psm_df.apply(lambda row: calculate_iRT_value(row.rt, iRT_m, iRT_t), axis=1)
psm_df['pep_seq'] = psm_df.apply(lambda row: create_unique_pep_seq(row), axis=1)

psm_df.sort_values("Score", inplace=True, ascending=False)

# filter out decoys
psm_df = psm_df[psm_df.isTT]

best_scores = []
for index, group in psm_df.groupby("pep_seq"):
    best_scores.append(group.head(1))

best_df = pd.concat(best_scores)
best_df2 = best_df.dropna(subset=['PSMID', 'PepSeq1', 'PepSeq2', 'LinkPos1', 'LinkPos2'])


lib_list = []
i = 1
for psm_index, psm in best_df2.iterrows():
    print("{}/{} Done.".format(i, best_df2.shape[0]))

    if not label_experiment:
        lbl = clName
    else:
        if check_for_isotope_labeling(psm.Crosslinker) == 0:
            lbl = clLightLabelName
        else:
            lbl = clHeavyLabelName

    entries = []
    xiAnn_json = get_annotation_xidb(psm)
    fragments = xiAnn_json['fragments']
    clusters = xiAnn_json['clusters']

    cl_residue_pair = [psm['ProteinLinkPos1'], psm['ProteinLinkPos2']]
    cl_residue_pair.sort()
    cl_residue_pair = [str(x) for x in cl_residue_pair]

    strippedSequence = strip_sequence(psm['PepSeq1'] + "_" + psm['PepSeq2'])
    labeledSequence = get_lbl_sequence(psm, lbl, label_both=True)

    if writeClSitesToModifiedSequence:
        modifiedSequence = get_lbl_sequence(psm, lbl, label_both=True)
    else:
        modifiedSequence = replace_mods(psm.PepSeq1 + "_" + psm.PepSeq2)

    entry_template = {
        "ProteinId": psm.Protein1 + "_" + psm.Protein2,
        "StrippedSequence": strippedSequence,
        "iRT": psm.iRT,
        "RT": psm.rt,
        "FragmentGroupId": psm.pep_seq,
        "PrecursorCharge": int(psm['match charge']),
        "PrecursorMz": psm['match mass']/psm['match charge'] + 1.00794,
        "ModifiedSequence": modifiedSequence,
        "IsotopeLabel": check_for_isotope_labeling(psm.Crosslinker),
        "scanID": psm['scan'],
        "run": psm.run,
        "searchID": psm.SearchID,
        #"peptide1": strip_sequence(psm['PepSeq1']),
        #"peptide2": strip_sequence(psm['PepSeq2']),
        "cl_residue_pair": "_".join(cl_residue_pair),
        "LabeledSequence": labeledSequence
    }

    for fragment in fragments:
        if not fragment['name'] == "P+P":
            for clusterId in fragment['clusterIds']:
                firstPeakId = clusters[clusterId]['firstPeakId']
                numbers = re.compile('\d+')
                entry = copy.deepcopy(entry_template)
                entry["FragmentCharge"] = clusters[clusterId]['charge']
                entry["FragmentType"] = fragment['name'][0]
                entry["FragmentNumber"] = numbers.findall(fragment['name'])[0]
                for cluster in fragment["clusterInfo"]:
                    if cluster['Clusterid'] == clusterId:
                        entry["FragmentMz"] = cluster['calcMZ']
                entry["RelativeFragmentIntensity"] = xiAnn_json['peaks'][firstPeakId]['intensity']
                # [x['intensity'] for x in xiAnn_json['peaks'] if clusterId in x['clusterIds']]
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

    entry_df = pd.DataFrame.from_records(entries)
    if len(entry_df) > 0:
        entry_df = entry_df[(entry_df.FragmentMz >= fragmentMzLowerLim) & (entry_df.FragmentMz <= fragmentMzUpperLim)]
    if len(entry_df) > 0:
        entry_df.drop_duplicates(inplace=True)
        entry_df.sort_values(by="RelativeFragmentIntensity", inplace=True, ascending=False)
        entry_df.RelativeFragmentIntensity = entry_df.apply(
                lambda row: row.RelativeFragmentIntensity / entry_df.RelativeFragmentIntensity.max(), axis=1)
        # ToDo: else error handling?
    lib_list.append(entry_df)
    i += 1
lib_df = pd.concat(lib_list)

# filter out neutral loss fragments
if not includeNeutralLossFragments:
    lib_df = lib_df[~lib_df.LossyFragment]

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


with open(output_path + "params.txt", "w") as text_file:
    psm_csv_path = baseDir + "psm_csv/"
    text_file.write("iRT_m : %s\n" % iRT_m)
    text_file.write("iRT_t : %s\n" % iRT_t)
    text_file.write("proteinId : %s\n" % proteinId)
    text_file.write("includeNeutralLossFragments: %s\n" % includeNeutralLossFragments)
    text_file.write("writeClSitesToModifiedSequence : %s\n" % writeClSitesToModifiedSequence)
    text_file.write("clName : %s\n" % clName)
    text_file.write("label_experiment : %s\n" % label_experiment)
    text_file.write("clLightLabelName : %s\n" % clLightLabelName)
    text_file.write("clHeavyLabelName : %s\n" % clHeavyLabelName)
    text_file.write("deuteriumCount : %s\n" % deuteriumCount)
    text_file.write("fragmentMzLowerLim : %s\n" % fragmentMzLowerLim)
    text_file.write("fragmentMzUpperLim : %s\n" % fragmentMzUpperLim)
    text_file.write("nClContainingFragments : %s\n" % nClContainingFragments)
    text_file.write("nLinearFragments : %s\n" % nLinearFragments)
    text_file.write("in files mzML:\n%s\n" % "\n".join(os.listdir(mzml_path)))
    text_file.write("in files PSM CSV:\n%s\n" % "\n".join(os.listdir(psm_csv_path)))
