Single-cell RNA-seq annotation with lineage barcodes
This repository contains code and annotation tables for associating 10x Genomics single-cell identifiers with lineage barcodes and clone classifications. The barcode-processing workflow generates counts for combinations of 10x cell barcodes and random DNA barcodes from sequencing reads. The anno1.txt–`anno4.txt` files associate single-cell data from four mice with clone output information.

Related preprint: https://www.biorxiv.org/content/10.1101/2025.02.23.639689v2


Part I — Repository contents
File	Description
10xBCコード.txt	Commands and Python snippets for trimming and combining FASTQ sequences and extracting/clustering barcodes with Bartender. This is a workflow note rather than a single executable script.
anno1.txt–`anno4.txt`	Cell-to-clone annotation tables for mice 1–4.
ex1_clone_census.py	Python script for assigning barcode families, filtering cross-clone barcode observations, identifying putative doublets, and producing a per-clone census.
HarmonyPCAanalysis.R	R analysis of clone-level distributions in Harmony embedding, neighbourhood enrichment, MolO comparisons, and differential expression.
README.md	Documentation.


In the annotation files, 10x denotes the single-cell identifier, BC the clone identifier, and Cluster the clone classification: c1 = P-HSCs; c2 = PRG-HSCs; c3 = PRGB-HSCs; c4 = balanced-HSCs.


Part II — System requirements
Software dependencies
Component	Required software / library	Version used or tested
FASTQ trimming	FASTX-Toolkit (fastx_trimmer)	To be confirmed for this analysis
Barcode extraction and clustering	Bartender (bartender_extractor_com, bartender_single_com)	To be confirmed for this analysis
FASTQ merging	Python (standard library)	3.13.5 (Anaconda base)
Clone census	Python	3.13.5 (Anaconda base)
Clone census	NumPy	To be confirmed for this analysis
Clone census	pandas	2.2.3 (Anaconda base)
Clone census	SciPy	1.15.3 (Anaconda base)
Harmony-based analysis	R	4.5.1
Harmony-based analysis	Seurat	5.5.0
Harmony-based analysis	Matrix	1.7.5
Harmony-based analysis	RANN	2.6.2
Harmony-based analysis	edgeR	Not installed in the checked R environment; version used for the original analysis is unconfirmed
Command-line execution	GNU Bash	3.2.57(1)-release


The specific Python dependencies of ex1_clone_census.py are numpy, pandas, scipy.cluster.hierarchy, and scipy.spatial.distance. There is no requirement to install Python's built-in sys module separately.
Operating systems and versions tested
The following operating system and software versions were used for this analysis:
- Operating system: macOS Sonoma 14.6.1 (build 23G93)
- Python: 3.13.5 (Anaconda base environment)
- pandas: 2.2.3 (Anaconda base environment)
- SciPy: 1.15.3 (Anaconda base environment)
- R: 4.5.1
- Seurat: 5.5.0
- Matrix: 1.7.5
- RANN: 2.6.2
- edgeR: Not installed in the checked R environment (original analysis version unconfirmed)
- GNU Bash: 3.2.57(1)-release
- Other tested operating systems: None documented
Versions of NumPy, FASTX-Toolkit, and Bartender remain to be documented. The checked R environment does not have edgeR installed; whether an alternate R environment was used for the edgeR-dependent analysis is not established. The fastx_env and treeenv environments used for the separate CPTracer pipeline are not listed here.
Hardware
Hardware requirements and typical execution time have not been documented. Memory requirements may depend on the number of cells and distinct barcodes; ex1_clone_census.py creates a dense cells-by-barcodes array and a barcodes-by-barcodes co-occurrence matrix.No additional installation is required for the code provided in this repository, apart from the third-party software listed under System requirements.


Part III — 10x barcode and lineage barcode processing
The following steps describe the commands in 10xBCコード.txt. Run them in a directory containing corresponding Read1.fastq and Read2.fastq files. Adapt file names to your data. Commands are presented in execution order; the two Python merging blocks in the original note must be run as Python code (they are not shell commands).
Step 1 — Extract read segments
fastx_trimmer -f 1 -l 16 -Q 33 -i Read1.fastq -o 10x-1.fastq
fastx_trimmer -f 37 -l 62 -Q 33 -i Read2.fastq -o BC.fastq
fastx_trimmer -f 32 -l 36 -Q 33 -i Read2.fastq -o ATTAA.fastq
These commands retain bases 1–16 from read 1, 37–62 from read 2, and 32–36 from read 2, respectively. The -Q 33 option specifies Phred+33 quality encoding.
Step 2 — Concatenate corresponding FASTQ sequences
The original workflow concatenates the ATTAA.fastq sequence/quality fields with those from 10x-1.fastq to make 10x-2.fastq, then concatenates 10x-2.fastq with BC.fastq to make 10xBC.fastq. The record identifiers from the first input are retained. Example equivalent Python implementation:
from itertools import zip_longest

def fastq_records(path):
    with open(path) as fh:
        while True:
            record = [fh.readline() for _ in range(4)]
            if not record[0]:
                return
            if any(not line for line in record):
                raise ValueError(f"Incomplete FASTQ record: {path}")
            yield [line.rstrip("\n\r") for line in record]

def concatenate_fastq(left_path, right_path, output_path):
    with open(output_path, "w") as out:
        for left, right in zip_longest(fastq_records(left_path), fastq_records(right_path)):
            if left is None or right is None:
                raise ValueError("FASTQ files have different numbers of records")
            if left[0].split()[0] != right[0].split()[0]:
                raise ValueError("FASTQ record identifiers do not match")
            if len(left[1]) != len(left[3]) or len(right[1]) != len(right[3]):
                raise ValueError("Sequence and quality lengths differ")
            out.write(f"{left[0]}\n{left[1]}{right[1]}\n{left[2]}\n{left[3]}{right[3]}\n")

concatenate_fastq("ATTAA.fastq", "10x-1.fastq", "10x-2.fastq")
concatenate_fastq("10x-2.fastq", "BC.fastq", "10xBC.fastq")
This is an explanatory rewrite of the merging logic rather than a separately validated program. The original source code is in 10xBCコード.txt.
Step 3 — Extract and cluster barcode sequences with Bartender
bartender_extractor_com -f 10xBC.fastq -o 10xBC -q '?' -p 'ATTAA[37]CCTAG' -m 2
bartender_single_com -f 10xBC_barcode.txt -o 10xBC -d 5
The barcode pattern and mismatch/distance parameters above reproduce the original workflow note. Check the generated files and ensure this pattern matches the relevant read design.


Part IV — Per-clone census (ex1_clone_census.py)
Step 1 — Prepare the input table
The standalone script expects a tab-delimited long-format file, with a header and four columns in this order:
Column	Meaning
cell	Single-cell identifier
BC	Barcode identifier
UMI	UMI count for the cell–barcode observation
RelativetoMax	UMI count divided by the maximum UMI count among barcodes in the same cell


The script validates RelativetoMax = UMI / max(UMI per cell) to within an absolute tolerance of 1e-6 for at least 99.9% of rows. The example input Ex1.txt is not included in the repository; it must be supplied separately. The script does not itself transform the Bartender output into this table.
Step 2 — Run the census script
python ex1_clone_census.py Ex1.txt Ex1_clone_census.tsv
The output is a tab-delimited per-clone census file. If arguments are omitted, the script defaults to Ex1.txt and Ex1_clone_census.tsv.
Step 3 — Interpret the analysis
The script:
1. Constructs a cell-by-barcode matrix of relative-to-maximum UMI counts.
2. Groups co-occurring barcodes using Jaccard distance with average-linkage hierarchical clustering (distance cutoff 0.6).
3. Merges nested barcode modules when the containment and relative cell-count criteria are met.
4. Iteratively assigns exclusive barcode ownership and removes observations attributed to other clones.
5. Reassigns cells with no retained barcode observations to the owner of their strongest barcode.
6. Flags candidate doublets when a foreign clone's barcode coverage ratio is at least 0.70, and excludes these cells from the per-clone census.
7. Reports clone size and barcode-related statistics, including barcode counts, mean/median barcodes per cell, relative-to-maximum UMI statistics, and the minimum barcode detection rate.
The thresholds are specified as constants in the script and were not re-optimized here.


Part V — Citation
Please cite the associated study when using these materials:



