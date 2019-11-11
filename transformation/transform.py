from collections import Counter
import HTSeq, os, json, subprocess, sqlite3, time
import numpy as np
import pandas as pd

def generate_seq_stats(seqfile, header, table=None, fastqfile=True):
    '''
    This function creates the JSON-files table.j, hist.j, edges.j, which are the basis for the sequence statistics table and graph visualized in the Sequence distribution-tab.
    If no table object is provided, headers are created and a table object is returned with two columns, headers and values. 
    If a table object is provided, the function will add a new column to the table

    table: existing table (for adding a column)
    seqfile: path to sequencefile (fasta/fastq)
    header: name of column
    fastqfile: the function assumes a fastq file. "False" will accept fasta
    '''
    if not table:
        table = {'Statistic': ['Count (#)','Length (bp)','Over 100 bp','Over 500 bp','Over 1000 bp','Over 5000 bp','Over 10000 bp','Largest (bp)','Smallest (bp)','Average length (bp)', 'Median (bp)', 'N50']}

    # Parse sequencefile
    if fastqfile:
       seqlengths = [len(s[0]) for s in HTSeq.FastqReader(seqfile, raw_iterator=True)]
    else:
       seqlengths = [len(s[0]) for s in HTSeq.FastaReader(seqfile, raw_iterator=True)]

    # Calculate statistcs
    table[header] = []
    table[header].append(len(seqlengths))
    table[header].append(sum(seqlengths))
    table[header].append(len([x for x in seqlengths if x>100]))
    table[header].append(len([x for x in seqlengths if x>500]))
    table[header].append(len([x for x in seqlengths if x>1000]))
    table[header].append(len([x for x in seqlengths if x>5000]))
    table[header].append(len([x for x in seqlengths if x>10000]))
    table[header].append(max(seqlengths))
    table[header].append(min(seqlengths))
    table[header].append(np.mean(seqlengths))
    table[header].append(calculate_n50(seqlengths))

    # Create historgram data
    hist, edges = np.histogram(seqlengths, density=False, bins=int(max(seqlengths)/10))
    return (table, hist.tolist(), edges.tolist())


def calculate_n50(numlist):
    '''
    This function takes a list of integers and calculates the n50 value for the given list. Returns float.
    Used in generate_seq_stats() to calculate n50 of any sequence set

    numlist: list of integers
    '''
    numlist.sort()
    newlist = []
    for x in numlist:
        newlist += [x]*x
    # take the mean of the two middle elements if there are an even number
    # of elements. otherwise, take the middle element
    if len(newlist) % 2 == 0:
        medianpos = int(len(newlist)/2)
        return float(newlist[medianpos] + newlist[medianpos-1]) /2
    else:
        medianpos = int(len(newlist)/2)
        return float(newlist[medianpos])

def parse_functional_mga(seqfile, data=None):
    """
    This function parses predicted gene sequences and returns a dict of dict data structure shaped for Bokeh table visualization

    data: a dict with already existing data (returned from any other related parse_functional_X function)
    seqfile: predicted genes (aminoacid sequences)
    """
    if not data:
        data = {}

    # HTseq implentation goes here
    for s in HTSeq.FastaReader(seqfile, raw_iterator=True):
        if not s[1] in data:
            data[s[1]] = {}
        data[s[1]]['seq'] = s[0]
    return data


def parse_functional_interpro(tsvfile, data=None):
    '''
    This function parses interpro tsv files and returns a dict of dict data structure shaped for Bokeh table visualization

    data: a dict with already existing data (returned from any other related parse_functional_X function)
    tsvfile: output tsv from interpro
    '''
    if not data:
        data = {}
    with open(tsvfile, "r") as tsv:
       for line in tsv:
           line = line.strip().split("\t")
           annotation = "{}|{}".format(line[4].strip(), line[5].strip())
           if not line[0].strip() in data:
               data[line[0].strip()] = {}
           data[line[0].strip()][line[3].strip()] = annotation

    return data

def parse_functional_to_sql(data, filename="csvdescriptions.sql"):
    '''
    This function writes an sqlite3 database from a transposed pandas dataframe for efficient searches in the visuliztion module

    data: pandas df from functional annotation parser
    filename: output filename
    '''
    df = pd.DataFrame(data).T
    conn = sqlite3.connect("output/csvdescriptions.sql")
    df.to_sql("functional", conn, if_exists='replace', index=True)

def parse_binning_maxbin(renamed, summary):
    '''
    This function parses the summary file from maxbin as well as information regarding fasta files. Returns a dict

    dirpath: path to the maxbin output directory
    summary: path to the summary file produced by the maxbin tool
    '''
    dfdata = pd.read_csv(summary, sep='\t')
    dfdata = dfdata.to_dict()
    
    dfdata['Bin annotation'] = {}
    taxonomy = ['taxonomy','aef','eaf','aeffea','feafea','feafea','feafea','feafe','fe']

    for count, key in enumerate(dfdata['Bin name'].keys(), 0):
        # Should be count once taxonomy parses is implemented
        dfdata['Bin annotation'][key] = taxonomy[0]
    return dfdata


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser("This script iterates through META-pipe result files and generates JSON data files needed for visualization. writes to data/")
    # Preprocessing
    parser.add_argument("--r1", help="Path to input R1 fastq")
    parser.add_argument("--r2", help="Path to input R2 fastq")
    parser.add_argument("--smerged", help="Path to input merged fastq")
    parser.add_argument("--sr1", help="Path to R1 remainder after merging (fastq)")
    parser.add_argument("--sr2", help="Path to R2 remainder after merging (fastq)")
    parser.add_argument("--str1", help="Path to QC processed R1")
    parser.add_argument("--str2", help="Path to QC processed R2")
    parser.add_argument("--stmerged", help="Path to QC processed merged fastq")
    parser.add_argument("--contigs", help="Path to contigs from assembly")
    # Taxonomic classification
    parser.add_argument("--rna16s", help="Path to 16s rRNA sequence file")
    parser.add_argument("--mapseq", help="Path to mapseq output file")
    parser.add_argument("--kaiju", help="Path to kaiju")
    # Functional annotation
    parser.add_argument("--mga", help="Path to predicted genes (amino acid seqs)")
    parser.add_argument("--interpro", help="Path to interpro annotation (tsv format)")
    parser.add_argument("--uniref", help="Path to uniref output")
    parser.add_argument("--mar", help="Path to MarRef output")
    # Binning
    parser.add_argument("--maxbinrenamed", help="Path to maxbin renaming conventions file")
    parser.add_argument("--maxbinsummary", help="Path to maxbin summary file")
    args = parser.parse_args()

    if not os.path.isdir("output"):
        os.mkdir("output")

    ### Sequence distribution
    if os.path.isfile(args.r1) and os.path.isfile(args.r2):
        seqhist = {}
        seqedge = {}
        table, seqhist['R1'], seqedge['R1'] = generate_seq_stats(args.r1, 'R1')
        table, seqhist['R2'], seqedge['R2'] = generate_seq_stats(args.r2, 'R2', table)
        table, seqhist['sMerged'], seqedge['sMerged'] = generate_seq_stats(args.sr1, 'sMerged', table)
        table, seqhist['sR1'], seqedge['sR1'] = generate_seq_stats(args.sr2, 'sR1', table)
        table, seqhist['sR2'], seqedge['sR2'] = generate_seq_stats(args.smerged, 'sR2', table)
        table, seqhist['stMerged'], seqedge['stMerged'] = generate_seq_stats(args.str1, 'stMerged', table)
        table, seqhist['stR1'], seqedge['stR1'] = generate_seq_stats(args.str2, 'stR1', table)
        table, seqhist['stR2'], seqedge['stR2'] = generate_seq_stats(args.stmerged, 'stR2', table)
        table, seqhist['Contigs'], seqedge['Contigs'] = generate_seq_stats(args.contigs, 'Contigs', table, False)
        
        with open('output/table.j', 'w') as output:
            json.dump(table, output)
        with open('output/hist.j', 'w') as output:
            json.dump(seqedge, output)
        with open('output/edges.j', 'w') as output:
            json.dump(seqhist, output)

    ### Taxonomic Classification
    # Both of these probably need conversions before Krona-scripts!
    if os.path.isfile(args.mapseq):
        pass
        #subprocess.run(['$HOME/Krona/KronaTools/scripts/ImportTaxonomy.pl', '-tax', '$HOME/taxonomy', '-o', '$HOME/output/mapseq.html', args.mapseq], shell=True)

    if os.path.isfile(args.kaiju):
        pass
        #subprocess.run(['$HOME/Krona/KronaTools/scripts/ImportTaxonomy.pl', '-tax', '$HOME/taxonomy', '-o', '$HOME/output/kaiju.html', args.kaiju], shell=True)

    ### Functional annotation
    data = None
    if os.path.isfile(args.mga):
        data = parse_functional_mga(args.mga)
    if os.path.isfile(args.mar):
        pass
    if os.path.isfile(args.interpro):
        data = parse_functional_interpro(args.interpro, data)
    if os.path.isfile(args.uniref):
        pass
    if data:
        parse_functional_to_sql(data)
    
    #### Binning
    if os.path.isfile(args.maxbinrenamed) and os.path.isfile(args.maxbinsummary):
        binning = parse_binning_maxbin(args.maxbinrenamed, args.maxbinsummary)
        with open('output/binning.j', 'w') as output:
            json.dump(binning, output)
        # Need (out).summary and (out).marker and fastafiles
