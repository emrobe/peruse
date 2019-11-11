import json, pprint
from collections import Counter
import HTSeq
import numpy as np
import pandas as pd
from scipy.stats.kde import gaussian_kde
from bokeh.sampledata.autompg import autompg as df

# jq -s . JSONFILE > JSONFILE to collect multiple objects into one array need for this parser
def get_dbs(data):
    tmp = []
    for entry in data:
        for feature in entry['features']:
            for annotation in feature['annotations']:
                if 'unirefId' in annotation:
                    tmp.append('UniRef50')
                if 'records' in annotation:
                    for record in annotation['records']:
                        if 'analysis' in record:
                            tmp.append(record['analysis'])
    dbs = list(set(tmp))
    return dbs

def calculate_n50(numlist):
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
        return newlist[medianpos]

def add_metrics(stats):
    metrics= {'Metric': [], 'Reads':[], 'Contigs':[], 'Genes': []}
    stats['genelengths'] = stats['genelengthsfull'] + stats['genelengthspartial']
    # Create a list with all genes...
    metrics['Metric'].append('Count (#)')
    metrics['Contigs'].append(len(stats['contiglengths']))
    metrics['Reads'].append(len(stats['readlengths']))
    metrics['Genes'].append(len(stats['genelengths']))
    # Total bp
    metrics['Metric'].append('Length (bp)')
    metrics['Contigs'].append(sum(stats['contiglengths']))
    metrics['Reads'].append(sum(stats['readlengths']))
    metrics['Genes'].append(sum(stats['genelengths']))
    # Over / under
    metrics['Metric'].append('Over 100 bp')
    metrics['Contigs'].append(len([x for x in stats['contiglengths'] if x>100]))
    metrics['Reads'].append(len([x for x in stats['readlengths'] if x>100]))
    metrics['Genes'].append(len([x for x in stats['genelengths'] if x>100]))
    metrics['Metric'].append('Over 500 bp')
    metrics['Contigs'].append(len([x for x in stats['contiglengths'] if x>500]))
    metrics['Reads'].append(len([x for x in stats['readlengths'] if x>500]))
    metrics['Genes'].append(len([x for x in stats['genelengths'] if x>500]))
    metrics['Metric'].append('Over 1000 bp')
    metrics['Contigs'].append(len([x for x in stats['contiglengths'] if x>1000]))
    metrics['Reads'].append(len([x for x in stats['readlengths'] if x>1000]))
    metrics['Genes'].append(len([x for x in stats['genelengths'] if x>1000]))
    metrics['Metric'].append('Over 5000 bp')
    metrics['Contigs'].append(len([x for x in stats['contiglengths'] if x>5000]))
    metrics['Reads'].append(len([x for x in stats['readlengths'] if x>5000]))
    metrics['Genes'].append(len([x for x in stats['genelengths'] if x>5000]))
    metrics['Metric'].append('Over 10000 bp')
    metrics['Contigs'].append(len([x for x in stats['contiglengths'] if x>10000]))
    metrics['Reads'].append(len([x for x in stats['readlengths'] if x>10000]))
    metrics['Genes'].append(len([x for x in stats['genelengths'] if x>10000]))
    # Max / Min
    metrics['Metric'].append('Largest (bp)')
    metrics['Contigs'].append(max(stats['contiglengths']))
    metrics['Reads'].append(max(stats['readlengths']))
    metrics['Genes'].append(max(stats['genelengths']))
    metrics['Metric'].append('Smallest (bp)')
    metrics['Contigs'].append(min(stats['contiglengths']))
    metrics['Reads'].append(min(stats['readlengths']))
    metrics['Genes'].append(min(stats['genelengths']))
    # Average, average-like
    metrics['Metric'].append('Average length (bp)')
    metrics['Contigs'].append(np.mean(stats['contiglengths']))
    metrics['Reads'].append(np.mean(stats['readlengths']))
    metrics['Genes'].append(np.mean(stats['genelengths']))
    metrics['Metric'].append('Median (bp)')
    metrics['Contigs'].append(np.median(stats['contiglengths']))
    metrics['Reads'].append(np.median(stats['readlengths']))
    metrics['Genes'].append(np.median(stats['genelengths']))
    metrics['Metric'].append('N50')
    metrics['Contigs'].append(calculate_n50(stats['contiglengths']))
    metrics['Reads'].append(calculate_n50(stats['readlengths']))
    metrics['Genes'].append(calculate_n50(stats['genelengths']))
    return metrics

# Hot to csv from dict ###
# data = {'row1': [1,2,3]}
# >>> pd.DataFrame.from_dict(data, orient='index', columns=['A','B','C'])

def count(data, sampleid):

    stats = {
            'sampleid': sampleid,
            'contigcount': len(data),
            'contiglengths': [],
            'genecount': 0,
            'genelengthsfull': [],
            'genelengthspartial': [],
            'interprocounts': [],
            'goterms':[],
            'dbs' : get_dbs(data),
            'csvacc': {},
            'csvdesc': {}
            }
    # Contigs
    for entry in data:
            stats['contiglengths'].append(len(entry['sequence']['value']))
            for feature in entry['features']:
                    # Genes
                    if 'location' in feature:
                            stats['genecount'] += 1
                    #pp.pprint(feature)
                    if feature['annotations'][0]:
                        unique_gene_name = str(entry['id'])+'_'+str(feature['annotations'][0]['geneId'])
                        stats['csvacc'][unique_gene_name] = {}
                        stats['csvdesc'][unique_gene_name] = {}
                    # Geneprediction info, UniRef db info, Uniref blast info and records, containing interpro for some reason
                    for annotation in feature['annotations']:
                            # Store amino acid sequence in csv
                            if 'value' in annotation:
                                stats['csvdesc'][unique_gene_name]['seq'] = annotation['value']
                                stats['csvacc'][unique_gene_name]['seq'] = annotation['value']
                            # Special case, get Uniref50 annotation to stats['csv']
                            if 'unirefId' in annotation:
                                db = annotation['unirefId'].split("_")[0]
                                stats['csvdesc'][unique_gene_name][db] = annotation['name'].strip().replace("\t", " ")
                                stats['csvacc'][unique_gene_name][db] = annotation['unirefId'].strip().replace("\t", " ")
                            # Count complete and partial genes, add lengths to lists in stats
                            if 'geneId' in annotation:
                                genelength = (annotation['endPos'] - annotation['startPos'])
                                if annotation['completePartial'] == '11':
                                    stats['genelengthsfull'].append(genelength)
                                else:
                                    stats['genelengthspartial'].append(genelength)
                            #pp.pprint(annotation)
                            #print type(annotation)
                            # For some reason, interpro has its own 'records' list
                            if 'records' in annotation:
                                    #pp.pprint (annotation['records'])
                                    #print type(annotation['records'])
                                    
                                    
                                    counts = []
                                    for blob in annotation['records']:
                                            # Special case for all dbs, as we might want to tailor specific parsing
                                            # If coils, parse "Coil, Length X"
                                            if blob['analysis'] == 'Coils':
                                                stats['csvdesc'][unique_gene_name][blob['analysis']] = "{}, Length: {}".format(blob['signatureAccession'], (blob['stopLocation']-blob['startLocation']))
                                            # If Gene3D, parse "interProDescription|signatureAccession"
                                            if blob['analysis'] == 'Gene3D':
                                                if 'interProAccession' in blob:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = blob['interProDescription']+'|'
                                                else:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = ""
                                                stats['csvdesc'][unique_gene_name][blob['analysis']] += blob['signatureAccession']
                                            # If Hamap, parse "interProDescription|signatureAccession"
                                            if blob['analysis'] == 'Hamap':
                                                if 'interProAccession' in blob:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = blob['interProDescription']+'|'
                                                else:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = ""
                                                stats['csvdesc'][unique_gene_name][blob['analysis']] += blob['signatureAccession']
                                            # If PIRSF, parse "interProDescription|signatureAccession"
                                            if blob['analysis'] == 'PIRSF':
                                                if 'interProAccession' in blob:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = blob['interProDescription']+'|'
                                                else:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = ""
                                                stats['csvdesc'][unique_gene_name][blob['analysis']] += blob['signatureAccession']
                                            # If PRINTS, parse "interProDescription|signatureAccession"
                                            if blob['analysis'] == 'PRINTS':
                                                if 'interProAccession' in blob:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = blob['interProDescription']+'|'
                                                else:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = ""
                                                stats['csvdesc'][unique_gene_name][blob['analysis']] += blob['signatureAccession']
                                            # If ProDom, parse "interProDescription|signatureAccession"
                                            if blob['analysis'] == 'ProDom':
                                                if 'interProAccession' in blob:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = blob['interProDescription']+'|'
                                                else:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = ""
                                                stats['csvdesc'][unique_gene_name][blob['analysis']] += blob['signatureAccession']
                                            # If ProSitePatterns, parse "interProDescription|signatureAccession"
                                            if blob['analysis'] == 'ProSitePatterns':
                                                if 'interProAccession' in blob:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = blob['interProDescription']+'|'
                                                else:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = ""
                                                stats['csvdesc'][unique_gene_name][blob['analysis']] += blob['signatureAccession']
                                            # If ProSiteProfiles, parse "interProDescription|signatureAccession"
                                            if blob['analysis'] == 'ProSiteProfiles':
                                                if 'interProAccession' in blob:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = blob['interProDescription']+'|'
                                                else:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = ""
                                                stats['csvdesc'][unique_gene_name][blob['analysis']] += blob['signatureAccession']
                                            # If SMART, parse "interProDescription|signatureAccession"
                                            if blob['analysis'] == 'SMART':
                                                if 'interProAccession' in blob:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = blob['interProDescription']+'|'
                                                else:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = ""
                                                stats['csvdesc'][unique_gene_name][blob['analysis']] += blob['signatureAccession']
                                            # If SUPERFAMILY, parse "interProDescription|signatureAccession"
                                            if blob['analysis'] == 'SUPERFAMILY':
                                                if 'interProAccession' in blob:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = blob['interProDescription']+'|'
                                                else:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = ""
                                                stats['csvdesc'][unique_gene_name][blob['analysis']] += blob['signatureAccession']
                                            # If TIGRFAM, parse "interProDescription|signatureAccession"
                                            if blob['analysis'] == 'TIGRFAM':
                                                if 'interProAccession' in blob:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = blob['interProDescription']+'|'
                                                else:
                                                    stats['csvdesc'][unique_gene_name][blob['analysis']] = ""
                                                stats['csvdesc'][unique_gene_name][blob['analysis']] += blob['signatureAccession']
                                            #if 'analysis' in blob:
                                            #    if 'interProAccession' in blob.keys():
                                            #        stats['csvacc'][unique_gene_name][blob['analysis']] = blob['interProAccession']+'|'
                                            #        stats['csvdesc'][unique_gene_name][blob['analysis']] = blob['interProDescription']+'|'
                                            #    else:
                                            #        stats['csvacc'][unique_gene_name][blob['analysis']] = ""
                                            #        stats['csvdesc'][unique_gene_name][blob['analysis']] = ""
                                            #    stats['csvacc'][unique_gene_name][blob['analysis']] += blob['signatureAccession']
                                            #    stats['csvdesc'][unique_gene_name][blob['analysis']] += blob['signatureDescription']
                                            # Gather all goterms for sample in stats['goterms']
                                            if blob['goAnnotations']:
                                                stats['goterms'] += blob['goAnnotations']
                                            #pp.pprint(blob['goAnnotations'])
                                            counts.append(blob['analysis'])
                                    counts = set(counts)
                                    counts = list(counts)
                                    #print counts
                                    stats['interprocounts'] += counts
    stats['interprocounts'] = Counter(stats['interprocounts'])
    stats['goterms'] = Counter(stats['goterms'])

    return stats

def count_reads(data):
    reads = {'readlengths': [], 'readcount': 0}
    reads['readlengths'] = [len(s[0]) for s in HTSeq.FastqReader(data, raw_iterator=True)]
    reads['readcount'] = len(reads['readlengths'])
    return reads

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser("This script iterates through META-pipe annotation files and returns a stats-object with visualization-ready statistics (schema not yet defined). The main function count() is used by other visualization scripts to gather data for plotting")
    parser.add_argument("--annotation", help="Annotations.j file from META-pipe. Rename file to <runID.j> to get correct run id in json object. Also, needs to be slurped using jq -s . first to collect all json objects in an array")
    parser.add_argument("--r1", help="R1 fastq")
    parser.add_argument("--r2", help="R2 fastq")
    args = parser.parse_args()

    with open(args.annotation, 'r') as data:
        data = json.load(data)
        sampleid = args.annotation.split("/")
        sampleid = sampleid[-1].split('.')
        sampleid = '.'.join(sampleid[:-1])
        stats = count(data, sampleid)

    with open(args.r1, 'r') as data:
        r1 = count_reads(data)

    with open(args.r2, 'r') as data:
        r2 = count_reads(data)
    
    r = {
         'readlengths': r1['readlengths'] + r2['readlengths'],
         'readcount': r1['readcount']+r2['readcount']
        } 
    # Add read data to stats
    stats.update(r)
    # Calculate extra metrics in stats
    metrics = add_metrics(stats)

    # Create export= {}, which is used for sequence distibution plots
    export = {}
    # Histogram data
    export['readhist'], export['readedges'] = np.histogram(stats['readlengths'], density=False, bins=int(max(stats['readlengths'])/10))
    export['contighist'], export['contigedges'] = np.histogram(stats['contiglengths'], density=False, bins=int(max(stats['contiglengths'])/10))
    export['genesfullhist'], export['genesfulledges'] = np.histogram(stats['genelengthsfull'], density=False, bins=int(max(stats['genelengthsfull'])/10))
 
    # Make sure array has data. Genelengthspartial might be empty due to config from META-pipe (MarCat)
    if stats['genelengthspartial']:
        export['genespartialhist'], export['genespartialedges'] = np.histogram(stats['genelengthspartial'], density=False, bins=int(max(stats['genelengthspartial'])/10))
    else:
        export['genespartialhist'], export['genespartialedges'] = np.array([0]),np.array([0,1])

    # Convert all ndarrays from numpy to regular lists because... why on earth is this not possible with python JSON omgomgomg
    for key in export.keys():
        export[key] = export[key].tolist()

    # KDE data
    export['readkdex'] = np.linspace(0,max(stats['readlengths']), 200).tolist()
    kde = gaussian_kde(stats['readlengths'])
    export['readkdey'] = kde(export['readkdex']).tolist()

    export['contigkdex'] = np.linspace(0,max(stats['contiglengths']), 200).tolist()
    kde = gaussian_kde(stats['contiglengths'])
    export['contigkdey'] = kde(export['contigkdex']).tolist()

    export['genekdex'] = np.linspace(0,max(stats['genelengthsfull']), 200).tolist()
    kde = gaussian_kde(stats['genelengthsfull'])
    export['genekdey'] = kde(export['genekdex']).tolist()


    export['sampleid'] = stats['sampleid']
    export['readstotal'] = len(stats['readlengths'])
    export['contigstotal'] = len(stats['contiglengths'])
    export['genesfullstotal'] = len(stats['genelengthsfull'])
    export['genespartialtotal'] = len(stats['genelengthspartial'])

    with open('data/metrics.j', 'w') as output:
        json.dump(metrics, output)
    with open('data/export.j', 'w') as output:
        json.dump(export, output)
    with open('data/csvaccessions.j', 'w') as output:
        json.dump(stats['csvacc'], output)
    with open('data/csvdescriptions.j', 'w') as output:
        json.dump(stats['csvdesc'], output)

#    dfacc = pd.DataFrame(stats['csvacc'])
#    dfacc = dfacc.T
#    dfacc.to_hdf('data/annotation_acc.hdf', 'dfacc', format='t')

#    dfdesc = pd.DataFrame(stats['csvdesc'])
#    dfdesc = dfdesc.T
#    dfdesc.to_hdf('data/annotation_desc.hdf5', 'dfdesc')
