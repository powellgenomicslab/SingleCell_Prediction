---
title: Brie functioning
author: Jose Alquicira Hernandez
---

Input files

- Genome reference (fasta file)
- Genome annotation (gtf fle)
- A BigWig file


# Get splicing events

`brie-event` uses an annotation file to obtain splicing events such as:

- Alternative 5' splice sites
- Alternative 3' splice sites
- Mutually exclusive exons
- Retained introns

It generates GFF3 files with splicing events information

`brie-event-filter` retrieves gold-quality splicing events only. See constraint criteria on [BRIE: splicing events](https://brie-rna.sourceforge.io/manual.html#splicing-events).

# Get sequence features

`brie-factor` gets a set of short sequence feature from the annotation file.

> `brie factor` needs `bigWigSummary` to work, it may be downloaded from [UCSC website](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary). See below.

```bash
cd /shares/common/users/j.alquicira/software/brie-master
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary
--2017-04-27 18:43:28--  http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary
Resolving hgdownload.cse.ucsc.edu... 128.114.119.163
Connecting to hgdownload.cse.ucsc.edu|128.114.119.163|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 4359384 (4.2M) [text/plain]
Saving to: `bigWigSummary'

100%[==================================================================================================================================>] 4,359,384   1.45M/s   in 2.9s    

2017-04-27 18:43:32 (1.45 MB/s) - `bigWigSummary' saved [4359384/4359384]
```


```bash
chmod +x /shares/common/users/j.alquicira/software/brie-master
export PATH="/shares/common/users/j.alquicira/software/brie-master:$PATH"
```

A `csv` file with sequence features is generated.

# BRIE isoform estimate

`brie` quantifies the fraction of exon inclusion level

# Differential splicing

`brie-diff` detects differential splicing between many cells pair-wisely, including just two cells, by calculating Bayes factor.

![Brie pipeline](brie_pipeline.png)