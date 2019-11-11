# 1 TODO: Run custom scripts to transform mapseq / kaiju to KronaTools-ready input. Store in data/

# 2 Run transform
python transform.py \
--r1 ENV/tiny_r1.fastq \
--r2 ENV/tiny_r2.fastq \
--smerged ENV/tiny_r1.fastq \
--sr1 ENV/tiny_r1.fastq \
--sr2 ENV/tiny_r2.fastq \
--str1 ENV/tiny_r1.fastq \
--str2 ENV/tiny_r2.fastq \
--stmerged ENV/tiny_r1.fastq \
--contigs ~/Programming_tutorial_bio3323/Project/salm.fasta \
--mga ENV/PROKKA_03222019.faa \
--interpro ENV/exporter.fas.interpro.out.all \
--uniref ENV/genes.faa \
--mar ENV/notafile \
--mapseq ENV/notafile \
--kaiju ENV/notafile \
--maxbinrenamed ENV/maxbin/renamed \
--maxbinsummary ENV/maxbin/PIN52_maxbin.summary
