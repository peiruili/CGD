# CGD
[![Build Status](https://travis-ci.org/joemccann/dillinger.svg?branch=master)](https://travis-ci.org/joemccann/dillinger)

CGD is python based on-target scoring method for CRISPRi, CRISPRa, Cas9 and Cas12a sequences.They are consortium of scores such as CGDi, CGDa, CGD9 and CGD12a, which are based on ENLOR(ELastic Net Logistic Regression) weights. 

### Requirements

  - Python 2.7 or above
  - Vienna RNA package
  - Numpy (version > = 1.12.4)
  - math
  - reverse complement (The file is provided as "get_sequence" in CGD)
  
### Package installation

For Vienna RNA package see the instructions in their official site https://www.tbi.univie.ac.at/RNA/ for python installation
The rest packages could be installed by using ```pip```

### Help
To make our source code user friendly we have provided help function which could be implemented in following way
```
python CGD.py -h help
usage: CGD.py [-h] [-a] [-b] [-c] [-d] [-e]

working with CGD

optional arguments:
  -h, --help  show this help message and exit
  -a        The option generates Comprehensive score using Comprehesive function
  -b        The option generates CRISPRi score using CGDi function
  -c        The option generates CRISPRa score using CGDa function
  -d        The option generates CRISPR-Cas9 score using CGD9 function
  -e        The option generates CRISPR-Cas12a score using CGD12a function
  -f        The option generates CRISPR-Cas9 (non-canonical) score using CGD9NG function
```    

### Use 
CGD python file is executed on fasta file
- The input file should be in fasta format 
- The length of fasta file should be betweem 100nt and 10000nt
- The scores are generated with range [0,1], where the gRNA sequences with score greater than 0.5 are considered efficient.
- The output file are shown in tabular format which compromises of ID, Start,End, Strand, 30 or 34nt sequence (based on Cas system)   and on-target score 
``` python
    python CGD.py -a input.fa (Comprehensive Score)
    
    Example file: input.fa 
    
    >XM_030244935
    CGGCGCGGAGTGCGCCGGCGCGTCGTCGGGGACGCCGGGTCCAGGATCTTGCTAGGGAACCAGTGTTGTC
    GCGTCGTCCCGCCCCCTCGGGGCTTTTGCTCCCGTTAACTGTCGGCGGGGCAGGCTCCGCAGCGCAGGGC
    GACATGCCGGTGCGCTTCAAGGGATCACGAAAGAACCAGGTTTTATTTCAAAAAGAAGAGTTCCCTACCA
    TGACCCTCAGATTTCAAAATACCTGGAGTGGAACGGAACCGTCAGAAAGAAGGATACGCTTGTCCCACCA
    GAACCCCAGGCCTTTGGAACGCCAAAGCCACAAGAGGCTGAGCAAGGAGAAGATGCCAATCAAGAAGCAG
    TTCTCTCACTAGAGGCCTCCAGGGTTCCCAAGAGAACTCGGTCTCATTCTGCGGACTCGAGAGCTGAAGG
    GGTTTCAGACACTGTGGAAAAGCACCAGGGTGTCACGAGAAGCCATGCGCCAGTTAGCGCGGATGTGGAG
    CTGAGACCTTCCAGCAAACAACCTCTCTCCCAGAGCATAGATCCCAGGTTGGATAGGCATCTTCGTAAGA
    AAGCTGGATTGGCCGTTGTTCCCACGAATAATGCCTTGAGAAATTCTGAATACCAAAGGCAGTTTGTTTG
```
- Output (CGD.txt)
```
ID	Start	End	Strand	Sequence	CGDi	CGDa	CGD9	CGDNGR	CGD12a
XM_030244935	1877	1907	-	TGCTGGACGGATGACAGTAAGGCGGGGCAC	0.86	0.39	1.0	0.0	0.0
XM_030244935	1877	1907	-	TGCTGGACGGATGACAGTAAGGCGGGGCAC	0.86	0.39	1.0	0.0	0.0
XM_030244935	1877	1907	-	TGCTGGACGGATGACAGTAAGGCGGGGCAC	0.86	0.39	1.0	0.0	0.0
XM_030244935	1510	1540	+	AAGGGAGGCAGGCTTCCTACACCGAGGCTG	0.92	0.39	0.99	0.0	0.0
XM_030244935	1867	1897	-	ATGACAGTAAGGCGGGGCACATGGAGGCAG	0.96	0.7	0.99	0.0	0.0
XM_030244935	1301	1331	-	CCCAAGCCAGCTTCCTCCTGACCGAGGCGG	0.89	0.53	0.99	0.0	0.0
XM_030244935	1301	1331	-	CCCAAGCCAGCTTCCTCCTGACCGAGGCGG	0.89	0.53	0.99	0.0	0.0
XM_030244935	479	509	+	CGGATGTGGAGCTGAGACCTTCCAGCAAAC	0	0	0	0.98	0
XM_030244935	264	294	+	TACGCTTGTCCCACCAGAACCCCAGGCCTT	0	0	0	0.81	0
XM_030244935	228	262	+	ATACCTGGAGTGGAACGGAACCGTCAGAAAGAAG	0	0	0	0	0.98
XM_030244935	969	1003	+	TTGGACCCGGGTGAAGGAGAACCTGTCAAACCAG	0	0	0	0	0.97
```
```
  python CGD.py -b input.fa  (Exclusive for CRISPRi)
```
- Output (CGDi.txt)
```
ID	Start	End	Stand	Sequence	CGDi	
XM_030244935	65	95	+	TTGTCGCGTCGTCCCGCCCCCTCGGGGCTT	0.99
XM_030244935	1363	1393	+	GAGGAGCCCAGGGCGGAGGAGGACGGGAGA	0.98
XM_030244935	1377	1407	+	GGAGGAGGACGGGAGAGAGGAGAGAGGACA	0.86
XM_030244935	271	301	-	CGTTCCAAAGGCCTGGGGTTCTGGTGGGAC	0.85
```
```
  python CGD.py -c input.fa  (Exclusive for CRISPRa)
```
- Output (CGDa.txt)
```
ID	Start	End	Strand	Sequence	CGDa
XM_030244935	1363	1393	+	GAGGAGCCCAGGGCGGAGGAGGACGGGAGA	0.83
XM_030244935	2	32	+	GCGCGGAGTGCGCCGGCGCGTCGTCGGGGA	0.83
XM_030244935	1370	1400	+	CCAGGGCGGAGGAGGACGGGAGAGAGGAGA	0.72
XM_030244935	1362	1392	+	GGAGGAGCCCAGGGCGGAGGAGGACGGGAG	0.71
```
```
  python CGD.py -d input.fa  (Exclusive for CRISPR-Cas9)
```
- Output (CGD9.txt)
```
ID	Start	End	Strand	Sequence	CGD9
XM_030244935	1877	1907	-	TGCTGGACGGATGACAGTAAGGCGGGGCAC	1.0
XM_030244935	1867	1897	-	ATGACAGTAAGGCGGGGCACATGGAGGCAG	0.99
XM_030244935	1926	1956	+	AGACCCTGAGTTTCAGCACAACATGGGAAA	0.95
XM_030244935	569	599	-	TCAAGGCATTATTCGTGGGAACAACGGCCA	0.94
```
```
  python CGD.py -d input.fa  (Exclusive for CRISPR-Cas12a)
```
- Output (CGD12a.txt)
```
ID	Start	End	Strand	Sequence	CGD12a
XM_030244935	1831	1865	+	AAAATTTTGGACCGTCAGCCCAGCACCCCTGGGC	1.0
XM_030244935	198	232	-	GTATTTTGAAATCTGAGGGTCATGGTAGGGAACT	0.91
XM_030244935	758	792	+	AGTATTTGAAAGGAAACAGCAGTCTGGAGATGCT	0.9
XM_030244935	654	688	+	AGTGTTTGCATCCAATCAGTTCCAAGGCAATACA	0.89
```
```
  python CGD.py -d input.fa  (Exclusive for CRISPR-Cas12a)
```
- Output (CGD9NG.txt)
```
ID  Start	End	Strand	Sequence	CGD9NG
XM_0302449	1378	1408	+	GAGGAGGACGGGAGAGAGGAGAGAGGACAG	1.0
XM_0302449	1341	1371	+	GAGCACGAAGGAAGACACCCAGGAGGAGCC	0.99
XM_0302449	784	814	-	GAGATGCTGACTCCAGTAAAGAAGGCAGAT	0.98
XM_0302449	834	864	-	AGACATGGCGTCGGAAGACTCAGACCAGTC	0.8
```

- All the output files contains score greater than 0.5 


License
----
MIT
