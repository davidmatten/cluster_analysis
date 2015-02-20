#!/usr/bin/python

import sys, os
import argparse
import pandas as pd
import scipy
import scipy.stats

def infiles_to_dict(infile, biomeasures):
    data = {}
    df = pd.io.parsers.read_csv(infile, sep="\t")
    infile_ids = list(df["seq_id"])

    df2 = pd.io.parsers.read_csv(biomeasures, sep="\t")
    bio_ids = list(df2["seq_id"])
    intersection = list(set(bio_ids) & set(infile_ids))
#    print "OMITTING: " + str(list( set(bio_ids) - set(infile_ids)))
#    print "OMITTING: " + str(list( set(infile_ids) - set(bio_ids)))

    for id in intersection:
        data[id] = {"seq":None, "vl":None, "cd4":None}
    for row in df.iterrows():
        id = row[1][0]
        if id in intersection:
            seq = row[1][1]
            data[id]["seq"] = seq
    for row in df2.iterrows():
        id = row[1][0]
        if id in intersection:
            vl = row[1][1]
            cd4 = row[1][2]
            data[id]["vl"] = vl
            data[id]["cd4"] = cd4
    return data

def get_cons_seq(infile):
    fh = open(infile)
    for line in fh:
        if line.split("\t")[0].lower().strip() == "cons_c":
            return line.split("\t")[1]
    fh.close()
    return ""

def get_counts(data, pos):
    d = {}
    for k in data.keys():
        residue = data[k]["seq"][pos]
        if residue in d.keys():
            d[residue] += 1
        else:
            d[residue] = 1
    return d

def meets_criteria(pos_d):
    if len(pos_d.keys()) > 1:
        l = []
        for k, v in pos_d.items():
            l.append(v)
        m = max(l)
        l.remove(m)
        for i in l:
            if 5<=i:
                return True
    return False

def make_test_list(data, pos, col):
    d = {}
    for k, v in data.items():
        residue = v["seq"][pos]
        if residue not in d.keys():
            d[residue] = [v[col]]
        else:
            d[residue].append(v[col])
    return d

def main(infile, biomeasures):
    outfile = "analysis_outfile.csv"
    data = infiles_to_dict(infile, biomeasures)
    cons_seq = get_cons_seq(infile)

    df = pd.io.parsers.read_csv(infile, sep="\t")
    df2 = pd.io.parsers.read_csv(biomeasures, sep="\t")
    df3 = pd.merge(df, df2, on="seq_id", how="outer")

    fh = open(outfile, "w")
    fh.write("alignment position, VL p-value (t-test), VL t-value (t-test), VL p-value (mann-whitney), VL u-value (mann-whitney), mean log VL (non-cons), mean log VL (cons), number (non cons), number (cons), CD4 p-value (t-test), CD4 t-value (t-test),CD4 p-value (mann-whitney), CD4 u-value (mann-whitney), mean CD4 (cons), mean CD4 (non-cons)\n")
    seq_len = len(data[data.keys()[0]]["seq"])
    for pos in range(seq_len):
        outline = []
        pos_d = get_counts(data, pos)
        if meets_criteria(pos_d):
            d = make_test_list(data, pos, "vl")

            dfx = df3[ df3["seq"].str[pos] == cons_seq[pos]]
            dfy = df3[ df3["seq"].str[pos] != cons_seq[pos]]

            cons_like_vl = list(dfx["VL"])
            cons_like_vl = [x for x in cons_like_vl if str(x) != 'nan']
            cons_unlike_vl = list(dfy["VL"])
            cons_unlike_vl = [x for x in cons_unlike_vl if str(x) != 'nan']
            
            (t, p) =  scipy.stats.ttest_ind(cons_like_vl, cons_unlike_vl)
            (mw_u, mw_p) = scipy.stats.mannwhitneyu(cons_like_vl, cons_unlike_vl)
            outline.append(pos)
            outline.append(p)
            outline.append(t)
            outline.append(mw_p)
            outline.append(mw_u)
            outline.append( sum(cons_unlike_vl)*1. / len(cons_unlike_vl)*1. )
            outline.append( sum(cons_like_vl)*1. / len(cons_like_vl)*1.)
            outline.append( len(cons_unlike_vl) )
            outline.append( len(cons_like_vl) )

            cons_like_cd4 = list(dfx["log_CD4"])
            cons_like_cd4 = [x for x in cons_like_cd4 if str(x) != 'nan']
            cons_unlike_cd4 = list(dfy["log_CD4"])
            cons_unlike_cd4 = [x for x in cons_unlike_cd4 if str(x) != 'nan']

            (t, p) =  scipy.stats.ttest_ind(cons_like_cd4, cons_unlike_cd4)
            (mw_u, mw_p) = scipy.stats.mannwhitneyu(cons_like_cd4, cons_unlike_cd4)
            outline.append(p)
            outline.append(t)
            outline.append(mw_p)
            outline.append(mw_u)
            outline.append( sum(cons_unlike_cd4)*1. / len(cons_unlike_cd4)*1. )
            outline.append( sum(cons_like_cd4)*1. / len(cons_like_cd4)*1.)
           
            fh.write(",".join([str(x) for x in outline])+"\n")

    fh.close()
 #   print "end."



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Performs clustering analysis on tab delimited aligned protein sequences matching some biological measurement. This was developed using cd4 and log viral load. Example usage: python cluster_analysis.py -seqs='/path/to/sequence/file.txt' -bio_measure='/path/to/biological/measures.txt' -outf='/path/to/out.out'""")
    parser.add_argument('-seqs', '--sequences', type=str,
            help='The path to the tab delimited aligned protein sequence files. eg: "/path/to/sequences/files.txt"', required=True)
    parser.add_argument('-bio', '--biomeasures', type=str,
            help="The path to the tab delimited file containing biological measurements of the associated protein sequences in the sequences file.", required=True)
    args = parser.parse_args()
    
    seqs = args.sequences
    biomeasures = args.biomeasures

    if (seqs is None) or (biomeasures is None) or (outfile is None):
        print "Please specify all required arguments. infile, bio measures, outfile"
        sys.exit()

    if not os.path.exists(seqs):
        print "Please supply a valid path to the sequence source files"
        sys.exit()
    if not os.path.exists(biomeasures):
        print "Please supply a valid path to the bio measures source files"
        sys.exit()

    main(seqs, biomeasures)

