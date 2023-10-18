# FICLE (Full Isoform Characterisation from (Targeted) Long-read Experiments)

<img align="left" src="https://github.com/SziKayLeung/FICLE/assets/33493350/c998696c-6820-4b36-950f-603cbe2b93e3" width="400" height="450">

Developed for analysing short-read RNA-Seq data, existing tools for assessing alternative splicing events fail to capture the connectivity and complexity of isoforms detected from long-read sequencing, particularly those generated from targeted profiling where a deep sequencing coverage is achieved.

FICLE (**F**ull **I**soform **C**haracterisation from **L**ong-read sequencing **E**xperiments) is a tool subsequently developed to comprehensively annotate and accurately document the full reportoire of isoforms identifed from targeted full-length long-read sequencing experiment. It accurately assesses the occurrence of alternative splicing events by comparing splice sites (exon) coordinates between long-read-derived transcripts and reference transcripts, further easing visualisation of isoforms on the UCSC genome browser. 

<br><br>

___________
## Installation
The latest FICLE release (18/10/2023) is version 1.1.3. See our [wiki](https://github.com/SziKayLeung/FICLE/wiki/Running-FICLE) for installation instructions.

For information about previous releases (> v1.1.2) and features introduced in them, see the [version history](https://github.com/SziKayLeung/FICLE/wiki/Version-history).
___________
## Documentation
Please see the Wiki page for more details:
1. [Introduction](https://github.com/SziKayLeung/FICLE/wiki/Introduction)
2. [How does FICLE work?](https://github.com/SziKayLeung/FICLE/wiki/How-does-FICLE-work%3F)
3. [Running FICLE](https://github.com/SziKayLeung/FICLE/wiki/Running-FICLE)
4. [Understanding the output of FICLE](https://github.com/SziKayLeung/FICLE/wiki/FICLE-output)

___________
## Vignette
Demo input data and output results can be found in the [vignette](https://github.com/SziKayLeung/FICLE/tree/master/vignette) folder. 

___________
## Citations
SK.Leung, A.R.Jeffries,...E.Hannon, J.Mill. **Long-read transcript sequencing identifies differential isoform expression in the entorhinal cortex in a transgenic model of tau pathology**.
[doi: https://doi.org/10.1101/2023.09.20.558220](https://www.biorxiv.org/content/10.1101/2023.09.20.558220v1) 
