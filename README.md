Insider
=======

Insider is a program that identifies genomic positions where several aligned reads have the same start or end positions. The tool can be used to identify genomic positions for viral integration sites, or to find on- and off-target sites in CRISPR-Cas9 digested genomic DNA. In both cases, the sites are detected with high sensitivity and specificity by rapid processing of read alignments.


## Insider -- Usage and parameters

![parameters explained](https://github.com/UppsalaGenomeCenter/insider/raw/master/wiki/insider_wiki_s.jpg)

### Usage:
    ./insider.pl [options] --help -i -o -ots -offset -minp -mins -minc --silent

### Parameters:
    -i      Insider input file (SAM format).

    -o      Output file for Insider results.

    -ots    Set this flag for analysis of CRISPR-Cas9 cleavage sites in 
            off-target sequencing (OTS) data.

    -offset Offset between left-end and right-end soft clipped alignments in
            host genome, used when searching for insertion sites (default=10).

    -minp   Minimum number of soft clipped reads starting at the exact same
            position in the insertion site peak. This number of reads is
            required for both the left-end and right-end alignments (default=3).

    -mins   Minimum number soft clipped reads starting at the exact same
            position in the insertion site peak, when taking both left-end
            and right-end alignments into account (default=10).

    -minc   Minimum number of soft clipped bases required for a read to be
            considered in the Insider analysis (default=20).

    --silent Do not print status to stdout.


### Insider usage example: Cas9 cleavage site detection
 
In this example, we will use Insider for detection of CRISPR-Cas9 off-target sites in the human genome. The data was generated using the Nano-OTS protocol presented in [this study](https://doi.org/10.1101/2020.02.09.940486). The nanopore reads were aligned to the human reference genome (GRCh38) using [Minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778).
 
First, download the Insider git repository using the following command:

   `git clone https://github.com/UppsalaGenomeCenter/insider.git`

In the repository you will find both the executable file and an example SAM file. 

To execute Insider, simply run the following command:

    `./insider.pl -i Nano-OTS.example.minimap2.chr22.sam -o out.txt -ots`

The example SAM file contains reads for human chromosome 22, resulting from applying the Nano-OTS protocol for three gRNAs (ATXN10, MMP14 and NEK1) to DNA from the human cell line HEK293. The insider software identifies Cas9 cleavage sites for these gRNAs. For this example, the following results are obtained:

chr22	11211032	11211032	5<br>
chr22	11215498	11215498	7<br>
chr22	16937226	16937226	20<br>
chr22	45689924	45689924	5<br>
chr22	45794862	45794862	51<br>
chr22	16937206	16937206	5<br>
chr22	47094773	47094773	8<br>

The three first result columns display the chromosome, start and end positions for the Cas9 cleavage site (GRCh38 coordinates). The fourth column displays the number of reads in support of a Cas9 cleavage site. The coordinate chr22:45794862 corresponds to the on-target site for ATXN10 and all other coordinates are potential off-target sites. 




Copyright (C) 2015-2020 Ignas Bunikis, Adam Ameur

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
