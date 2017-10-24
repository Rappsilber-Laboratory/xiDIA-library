#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pandas as pd
import pymzml
import numpy as np
import json
import urllib.request
import datetime
import sys
import copy
import re
import os
import glob
import codecs


# ----Paths---

baseDir = "D:/Daten/Edinburgh/Fraenze/DIA paper 8 protein mix/Library/"
psm_csv_path = baseDir + "psm_csv/"
mzml_path = baseDir + "mzml/"
output_path = baseDir + "out/"

# iRT equation parameters y = m * x + t
# e.g. y= 0.6452 x + 47.469
iRT_m = 0.6452
iRT_t = 47.469

proteinId = "7PMix"
includeNeutralLossFragments = False

label_experiment = False
deuteriumCount = 4

nClContainingFragments = 15
nLinearFragments = 5

# ------------
# checks
if nClContainingFragments + nLinearFragments > 20:
    raise Exception("Total number of fragments should not exceed 20!")

# create psm_dataframe and merge xiIDs and gradient information
df_list = []
for csv in glob.glob(psm_csv_path + "*.csv"):
    temp_df = pd.read_csv(csv)
    df_list.append(temp_df)
psm_df = pd.concat(df_list, ignore_index=True)


# load mzmls
mzmls = []
for filename in os.listdir(mzml_path):
    if "mzML" in filename:
        mzmls.append(filename)

msruns = {}
for mzml in mzmls:
    msrun = pymzml.run.Reader(mzml_path + mzml)
    msruns[mzml[:-5]] = msrun


def get_rt_from_scanID(id1, run_name, msruns):
    msrun = msruns[run_name]
    # print run_name
    # print id1
    try:
        # rt = msrun[id1]['scan start time']
        rt = msrun[int(id1)]['scan start time']
    except:
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

    return '{}_{}{}-{}:{}'.format(pep1, pep2, linkpos1, linkpos2, charge)


def get_annotations(psm_row):
    """
    Runs the Xi annotator to get the fragment annotation.
    """
    if psm_row.SearchID > 10000:
        url = 'http://xi3.bio.ed.ac.uk/xiAnnotator/annotate/%s/85160-94827-76653-69142/%s/?peptide=%s&peptide=%s&link=%s&link=%s' % (
            int(psm_row.SearchID), int(psm_row.PSMID), psm_row.PepSeq1, psm_row.PepSeq2, int(psm_row.LinkPos1),
            int(psm_row.LinkPos2))
    else:
        url = 'http://129.215.14.63/xiAnnotator/annotate/%s/85160-94827-76653-69142/%s/?peptide=%s&peptide=%s&link=%s&link=%s' % (
            int(psm_row.SearchID), int(psm_row.PSMID), psm_row.PepSeq1, psm_row.PepSeq2, int(psm_row.LinkPos1),
            int(psm_row.LinkPos2))

    reader = codecs.getreader("utf-8")
    data = json.load(reader(urllib.request.urlopen(url)))
    return data


# def replace_HL_mods(sequence):
#     """ Removes heavy/light modifications of the cross-linker, e.g. bs3nh2 """
#     mods = ["bs3nh2", "bs3oh", "bs3-d4nh2", "bs3-d4oh"]
#     mods = "|".join(mods)
#     return replace_mods(re.sub('({})'.format(mods), '', sequence))


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
#     expmz = []
#     scanmz = []
#     for ii, row in psm_df.iterrows():
#
#         try:
#             scanmz.append(float(msruns[row["run"]][row["scan"]]["precursors"][0]["mz"]))
#             expmz.append(row["exp m/z"])
#         except:
#             pass
#
#     expmz = np.array(expmz)
#     scanmz = np.array(scanmz)
#     diff = scanmz - expmz
#     plt.hist(diff)


def invert_modifications(modstring, lth=True):
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
            modstring = modstring.replace(i, j)
    else:
        # replace light with heavy
        for i, j in zip(light_mods, heavy_mods):
            modstring = modstring.replace(j, i)

    return modstring


def calcIsotopeFragMz(mz, charge, cl_containing, label, label_mz_difference):
    """
    Adjusts the mass computation for the un-indentified labeled sequence
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
    return (unmatches)


def get_lbl_sequence(psm, lbl, label_both=True):
    """
    Creates the 'LabeledSequene' column for spectronaut.


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
    return (pepseq[:clpos + 1] + "[" + lbl + "]" + pepseq[clpos + 1:])


def check_AA(cl_pos1, modseq1, linkable, entry_template):
    """
    checks if the cl amino acid is a vlid one
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


# %%
# annotate psms with retention time
psm_df['rt'] = psm_df.apply(lambda row: get_rt_from_scanID(row["scan"], row["run"], msruns), axis=1)

# transform rts to iRT -> what about the nan columns?
psm_df['iRT'] = psm_df.apply(lambda row: calculate_iRT_value(row.rt, iRT_m, iRT_t), axis=1)
psm_df['pep_seq'] = psm_df.apply(lambda row: create_unique_pep_seq(row), axis=1)

psm_df.sort_values("Score", inplace=True, ascending=False)

best_scores = []
for index, group in psm_df.groupby("pep_seq"):
    best_scores.append(group.head(1))

best_df = pd.concat(best_scores)
best_df2 = best_df.dropna(subset=['PSMID', 'PepSeq1', 'PepSeq2', 'LinkPos1', 'LinkPos2'])

# %%
lib_list = []
i = 1
for psm_index, psm in best_df2.iterrows():
    #    if ("ox" in replace_mods(psm.PepSeq1 + "_" + psm.PepSeq2) )and ("bs3" in replace_mods(psm.PepSeq1 + "_" + psm.PepSeq2)):
    #        break
    if i % 50 == 0:
        print ("{}/{} Done.".format(i, best_df2.shape[0]))

    if check_for_isotope_labeling(psm.Crosslinker) == 0:
        lbl = "bs3-light"
    else:
        lbl = "bs3-heavy"

    entries = []
    xiAnn_json = get_annotations(psm)
    fragments = xiAnn_json['fragments']
    clusters = xiAnn_json['clusters']

    cl_residue_pair = [psm['ProteinLinkPos1'], psm['ProteinLinkPos2']]
    cl_residue_pair.sort()
    cl_residue_pair = [str(x) for x in cl_residue_pair]

    test = psm.PepSeq1 + "_" + psm.PepSeq2

    entry_template = {
        "ProteinId": psm.Protein1 + " _" + psm.Protein2,
        "StrippedSequence": strip_sequence(psm['PepSeq1'] + "_" + psm['PepSeq2']),
        "iRT": psm.iRT,
        "RT": psm.rt,
        "FragmentGroupId": psm.pep_seq,
        "PrecursorCharge": int(psm['match charge']),
        "PrecursorMz": psm['match mass']/psm['match charge']+1.00794,
        #"ModifiedSequence": replace_HL_mods(psm.PepSeq1 + "_" + psm.PepSeq2),
        "ModifiedSequence": replace_mods(psm.PepSeq1 + "_" + psm.PepSeq2),
        "IsotopeLabel": check_for_isotope_labeling(psm.Crosslinker),
        "scanID": psm['scan'],
        "run": psm.run,
        "searchID": psm.SearchID,
        #"peptide1": strip_sequence(psm['PepSeq1']),
        #"peptide2": strip_sequence(psm['PepSeq2']),
        "cl_residue_pair": "_".join(cl_residue_pair),
        "LabeledSequence": get_lbl_sequence(psm, lbl)
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
                entry["RelativeFragmentIntensity"] = xiAnn_json['peaks'][firstPeakId][
                    'intensity']  # [x['intensity'] for x in xiAnn_json['peaks'] if clusterId in x['clusterIds']]
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
        entry_df = entry_df[(entry_df.FragmentMz > 100) & (entry_df.FragmentMz < 1300)]
        entry_df.sort_values(by="RelativeFragmentIntensity", inplace=True, ascending=False)
        entry_df.RelativeFragmentIntensity = entry_df.apply(
                lambda row: row.RelativeFragmentIntensity / entry_df.RelativeFragmentIntensity.max(), axis=1)
        #else einfügen print i für error
    lib_list.append(entry_df)
    i += 1
lib_df = pd.concat(lib_list)

# filter out neutral loss fragments
if not includeNeutralLossFragments:
    lib_df = lib_df[~lib_df.LossyFragment]

filtered_fragments_list = []
min_CLfragments = 2
for index, group in lib_df.groupby(['FragmentGroupId', 'IsotopeLabel']):
    CLContaining = group[group.CLContainingFragment]
    nonCLContaining = group[~group.CLContainingFragment]
    topGroup = pd.concat([CLContaining.head(nClContainingFragments), nonCLContaining.head(nLinearFragments)])

    filtered_fragments_list.append(topGroup)


lib_df = pd.concat(filtered_fragments_list)

if label_experiment:
    missing_peptides = []
    print ("adding missing peptides...")

    deuterium_mass = 2.014102
    hydrogen_mass = 1.00794
    mzdiff = deuteriumCount * (deuterium_mass - hydrogen_mass)

    for i, group in lib_df.groupby('FragmentGroupId'):

        # if only heavy or light exists create the missing peptide
        if len(group.IsotopeLabel.unique()) < 2:

            missing_pep = copy.deepcopy(group)
            missing_pep['FragmentMz'] = missing_pep.apply(lambda row:
                calcIsotopeFragMz(row.FragmentMz, row.FragmentCharge, row.CLContainingFragment, row.IsotopeLabel, mzdiff), axis=1)

            pcharge = missing_pep["PrecursorCharge"].iloc[0]

            if missing_pep['IsotopeLabel'].mean() == deuteriumCount:
                missing_pep['IsotopeLabel'] = 0
                # count the modifications
                nmods = missing_pep['LabeledSequence'].iloc[0].count("bs3")
                missing_pep['PrecursorMz'] = missing_pep['PrecursorMz'].mean() - (mzdiff * nmods) / pcharge
                missing_pep['LabeledSequence'] = [invert_modifications(i.replace("[bs3-heavy]", "[bs3-light]"), lth=False)
                                                  for i in missing_pep['LabeledSequence']]

                # changes also the modifications strings...
                #missing_pep['ModifiedSequence'] = [invert_modifications(i, lth=False) for i in
                #                                   missing_pep['ModifiedSequence']]
                missing_pep['Artificial'] = "FromHeavy"

            else:
                nmods = missing_pep['LabeledSequence'].iloc[0].count("bs3")
                # generate the heavy version
                missing_pep['IsotopeLabel'] = 4
                missing_pep['PrecursorMz'] = missing_pep['PrecursorMz'].mean() + (mzdiff * nmods) / pcharge
                missing_pep['LabeledSequence'] = [invert_modifications(i.replace("[bs3-light]", "[bs3-heavy]"), lth=True)
                                                  for i in missing_pep['LabeledSequence']]
                # changes also the modifications strings...
                #missing_pep['ModifiedSequence'] = [invert_modifications(i, lth=True) for i in
                #                                   missing_pep['ModifiedSequence']]
                missing_pep['Artificial'] = "FromLight"

            # append any missing peptides to add them later to the df
            missing_peptides.append(missing_pep)

    lib_df = pd.concat(missing_peptides + [lib_df])

    # lib_df = complete_lib_df[["run", "scanID", "searchID", "ProteinId", "FragmentGroupId", "RT", "iRT",
    #                                    "Artificial", "IsotopeLabel", "peptide1", "peptide2",
    #                                    "LabeledSequence", "ModifiedSequence", "StrippedSequence",
    #                                    "CLContainingFragment", "FragmentCharge", "FragmentLossType",
    #                                    "FragmentMz", "FragmentNumber", "FragmentType",
    #                                    "LossyFragment", "PrecursorCharge", "PrecursorMz", "RelativeFragmentIntensity"]]

    # print ("copy to clipboard...")
    # complete_lib_df.to_clipboard(sep="\t")
    # print (get_all_mods(complete_lib_df["ModifiedSequence"]))


d = datetime.datetime.now()
d = datetime.date.strftime(d, "%m%d%y")
idstr = "_spectronaut_lbls"

# save the library
print ("Writing the library...")
lib_df.to_csv(output_path + proteinId + "_" + "_".join(
        [str(x) for x in lib_df.searchID.unique()]) + "nCL{}_nLin{}_lossy{}_lib{}.csv".format(
        nClContainingFragments, nLinearFragments, includeNeutralLossFragments, idstr), index=False)

#
# skyline_df = complete_lib_df[["FragmentMz", "StrippedSequence",
#                               "PrecursorMz", "ProteinId",
#                               "RelativeFragmentIntensity", "iRT"]]
#
# # change to skyline naming
# skyline_df.columns = ["ProductMz", "Peptide", "PrecursorMz", "Protein",
#                       "LibraryIntensity", "Tr_recalibrated"]
# skyline_df.to_csv("{}SkylineLib_{}.csv".format(output_path, complete_lib_df["ProteinId"].values[0]), index=False)
# complete_lib_df.to_csv("{}NormalLib_{}_{}.csv".format(output_path,
#                                                       complete_lib_df["ProteinId"].values[0],
#                                                       d), index=False)
#
# complete_lib_df.to_clipboard(sep="\t")

#
# # to_openms(best_df2, output_path)
# # ==============================================================================
# # extras tuff
# # ==============================================================================
# def to_openms(best_df2):
#     """
#     """
#     import pyopenms as oms
#     prot_ids = oms.ProteinIdentification()
#     peps = []
#
#     for ii, cl_pep in best_df2.iterrows():
#         prot_hit = oms.ProteinHit()
#         prot_hit.setAccession(cl_pep["Protein1"] + "_" + cl_pep["Protein2"])
#         prot_hit.setSequence((strip_sequence(cl_pep["PepSeq1"]) + "WW" +
#                               strip_sequence(cl_pep["PepSeq2"])))
#         prot_ids.insertHit(prot_hit)
#
#         # peptid ehit
#         pep_id = oms.PeptideIdentification()
#         pep_hit = oms.PeptideHit()
#         # pep_hit.setProteinAccessions([cl_pep.peptide1.name +
#         # "_CL_" + cl_pep.peptide2.name])
#         seq = oms.AASequence()
#         seq = seq.fromString(strip_sequence(cl_pep["PepSeq1"]) + "WW" + strip_sequence(cl_pep["PepSeq2"]), 0)
#         pep_hit.setSequence(seq)
#         pep_hit.setScore(cl_pep["Score"])
#         pep_hit.setCharge(int(cl_pep["exp charge"]))
#         # meta values
#         pep_hit.setMetaValue("MZ", cl_pep["exp m/z"])
#         pep_hit.setMetaValue("RT", cl_pep["rt"] * 60)
#         pep_hit.setMetaValue("seq1", cl_pep["PepSeq1"])
#         pep_hit.setMetaValue("seq2", cl_pep["PepSeq2"])
#         pep_hit.setMetaValue("clsite1", cl_pep["LinkPos1"])
#         pep_hit.setMetaValue("clsite2", cl_pep["LinkPos2"])
#         pep_hit.setMetaValue("scanID", cl_pep.scan)
#         pep_hit.setMetaValue("xiscore", cl_pep.Score)
#
#         # pep_id.insertHit(pep_hit)
#         pep_id.setRT(cl_pep["rt"] * 60)
#         pep_id.setMZ(cl_pep["exp m/z"])
#         pep_id.insertHit(pep_hit)
#         peps.append(pep_id)
#
#     id_file = oms.IdXMLFile()
#     # process
#     id_file.store(output_path + "identifications_OMS.idxml", [prot_ids], peps)
#
#     # create a featureXML
#     # iteate over map
#     featurexml = oms.FeatureXMLFile()
#     seeds = oms.FeatureMap()
#
#     for ii, cl_pep in best_df2.iterrows():
#         myfeat = oms.Feature()
#         myfeat.setMZ(cl_pep["exp m/z"])
#         myfeat.setRT(cl_pep["rt"] * 60)
#         myfeat.setMetaValue("seq1", cl_pep.PepSeq1)
#         myfeat.setMetaValue("seq2", cl_pep.PepSeq2)
#         myfeat.setMetaValue("label", cl_pep.PepSeq1 + "---" +
#                             cl_pep.PepSeq2)
#         myfeat.setPeptideIdentifications([pep_id])
#         myfeat.setIntensity(10 ** 5)
#         myfeat.setCharge(int(cl_pep["exp charge"]))
#         seeds.push_back(myfeat)
#
#     featurexml.store(output_path + "identifications_OMS.featurexml", seeds)
#
#
#     # fragment group - link: '{}_{}{}-{}:{}'.format(pep1, pep2, linkpos1, linkpos2, charge)
#
#     # # filter out precursors out of mz range
#     # upper_precursor_mz = 1200
#     # lower_precursor_mz = 500
#     # limited_df = lib_df[(lib_df.PrecursorMz > lower_precursor_mz) & (lib_df.PrecursorMz < upper_precursor_mz)]
#     #
#     # limited_df.to_csv("libs/C3_lib.csv", index=False)
#     # limited_df.to_csv("libs/"+"_".join(limited_df.searchID.unique())+"_lib.csv", index=False)
#     #
#     #
#     #
#     #
#     # # --test stuff-----
#     # lib_old = pd.read_csv("libs/HSA_BS3_injection_replica_-DMSO_01.08.16_PSM_xiFDR1.0.6.16_lib.csv", index_col=False)
#     # lib_new = pd.read_csv("libs/4724_lib.csv", index_col=False)
#     #
#     # lib_old.head(5)
#     # lib_new.head(5)
#     #
#     # lib_combined = pd.read_csv("libs/4884_4724_lib.csv", index_col=False)
#     #
#     # lib_single = pd.read_csv("libs/HSA_BS3_injection_replica_-DMSO_01.08.16_PSM_xiFDR1.0.6.16_lib.csv", index_col=False)
#     # lib_single["searchID"] = "4724"
#     # lib_single.FragmentGroupId = lib_single.apply(lambda row: row.FragmentGroupId[:-1], axis=1)
#     #
#     # lib_combined[lib_combined.searchID == 4724]
#     # overlap = lib_combined[lib_combined.FragmentGroupId.isin(lib_single.FragmentGroupId)]
#     #
#     # #overlap["test"] = overlap.apply(lambda row: lib_single[lib_single.FragmentGroupId == row.FragmentGroupId].FragmentMz, axis=1)
#     #
#     # overlap["ID"] = overlap.apply(lambda row: "_".join([str(row["FragmentGroupId"]), str(row["FragmentNumber"]), str(row["FragmentType"])]), axis=1)
#     # lib_single["ID"] = lib_single.apply(lambda row: "_".join([str(row["FragmentGroupId"]), str(row["FragmentNumber"]), str(row["FragmentType"])]), axis=1)
#     #
#     # len(set(overlap["ID"]) - set(lib_single["ID"]))
#     #
#     # joined = overlap.join(lib_single[["FragmentMz", "ID"]], on="ID", rsuffix="_single")
#     #
#     # joined[["FragmentMz_single", "FragmentMz"]]
#     # #
#     # # import seaborn as sns
#     # # import matplotlib.pyplot as plt
#     # # fig, ax = plt.subplots()
#     # #
#     # # sns.distplot(lib_df.FragmentMz, bins=20, kde=False, ax=ax)
#     # # ax.set_xticks(np.arange(100, 2100, 100))
#     # # ax.set_title(psm_csv[6:-4]+"_transitions")
#     # #
#     # #
#     # # fig, ax = plt.subplots()
#     # # sns.distplot(lib_df.PrecursorMz, bins=20, kde=False, ax=ax)
#     # # ax.set_xticks(np.arange(400, 1600, 100))
#     # # ax.set_title(psm_csv[6:-4]+"_precursor")
#     # #