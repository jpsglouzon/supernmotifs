# The super-n-motifs model #

**Comparing RNA secondary structures** of arbitrary size uncovers **structural patterns** that can provide a better understanding of **RNA functions**. However, performing fast and accurate **secondary structure comparisons** is **challenging** when we take into account the **RNA configuration (i.e., linear or circular)**, the presence of **pseudoknot** and **G-quadruplex (G4)** motifs and the **increasing number of secondary structures** generated by high-throughput probing techniques. To address this challenge, we propose the **super-n-motifs model**, based on a latent analysis of enhanced motifs comprising not only basic motifs but also adjacency relations. The super-n-motifs model computes a vector representation of secondary structures as linear combinations of these motifs.

* Download **[here](https://github.com/jpsglouzon/supernmotifs/releases)** the ubuntu 64bit executable and follow the 'how to use it' instructions to compare secondary structures. 

* Or follow the instructions below to compile and use the Super-n-motifs program on almost every computer platform.

## How to compile the Super-n-motifs program ? ##

* Download the source code **[here](https://github.com/jpsglouzon/supernmotifs/zipball/master)** and unzip.

* Open the terminal, `cd path_to_supernmotifs_program` to access the super-n-motifs program folder then compile it by running the command `make`.

* The executable file named 'supernmotifs' can be found in `path_to_supernmotifs_program`.

## How to use it? ##

* **Compare secondary structures** by calling: 
```
/path_to_supernmotifs_program/supernmotifs -i fileInDb -o folderOfResults
```

* The Super-n-motifs program takes as input a file of **rna secondary structures** in **dot-bracket** format (-i fileInDb):
```
>RNA1
GCCCCGCUGAUGAGGUCAGGGAAAACCGAAAGUGUCGACUCUACGGGGC
((((((.......((((......))))...((((....)))).))))))
```
It ouputs a **dissimilarity matrix** and various stats by default (-o folderOfResults) .
For further options: `./path_to_supernmotifs_program/supernmotifs -h`

* The super-n-motifs model supports circular RNA, pseudoknots and g-quadruplexes. **Circular RNAs** require to add `c_` to the header of each RNA :  `>c_RNA1`. Base pairs involved in **pseudoknots** are typically represented by special characters `{}`, `<>`, `[]`, and alphabets such as `Aa` or `Bb` : `..AA..aa`. Interacting guanines in **g-quadruplexes** are represented by `+` : `..((..++.++.++..++.)).`.

## How to cite the Super-n-motifs model ? ##

Glouzon JS, Perreault JP, Wang S. The super-n-motifs model: a novel
alignment-free approach for representing and comparing RNA secondary structures. Bioinformatics. 2017 Jan 14. pii: btw773. doi: 10.1093/bioinformatics/btw773

## Licence ##

The Super-n-motifs model is released under the terms of the GNU GPL licence. For further informations, please see the LICENCE file of the repository.


