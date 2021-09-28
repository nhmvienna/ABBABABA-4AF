# README

## The ABBA-BABA test

The original ABBA-BABA test for gene flow, conceived by Nick Patterson (Green et al. 2011), is based on a hypothetical phylogeny with one distant outgroup (H4) and three ingroup samples (H1,H2,H3,H4).

          /\
         /  \
        /\   \
       /  \   \
      /\   \   \
     /  \   \   \
    H1  H2  H3  H4

    A   B   B   A
    B   A   B   A

This approach uses haplotype sequencing data and counts the occurence of either of two genotype patterns "A-B-B-A" and "B-A-B-A" in the four samples, which are not conistent with the tree topology (see below) and may be indicative of incomplete lineage sorting and/or gene flow. In the case of incomplete lineage sorting only, both ABBA and BABA should occur at equal proportions.

<img src=/images/DG.png height="55" />

Thus count(ABBA)-count(BABA), scaled by the total number of counts (_D_<sub>G</sub>; see above) should equal 0. However, if there is geneflow from H3 to H2, we expect an excess of ABBA, resulting in positive values of _D_<sub>G</sub>.

Equivalent to Green _et al._’s (2010) test, Durand _et al._ (2011; _D_<sub>D</sub>) and Soraggi _et al._ (2018; _D_<sub>S</sub>) recently proposed using allele frequencies instead of counts. Here, _p_<sub>ij</sub> is the frequency of the derived allele B at the _i_<sup>th</sup> locus in the _j_<sup>th</sup> population in the generic phylogeny (P1, P2, P3 or P4).

<img src=/images/DD.png height="50" />

<img src=/images/DS.png height="50" />

## ABBABABA-4AF

The python script "ABBABABA-4AF.py" is a new implementation of the ABBA-BABA test for gene flow based on the latter two approaches and calculates both genome-wide and window-wise estimates of _D_ based on a given (fixed) number of SNPs per window. A 2D matrix of allele frequencies is used as input where columns represent samples and rows SNP positions. Following Green _et al._ (2011), the software also calculates _z_-scores based on jackknifing using a blocked, even _m_-delete jackknife procedure, where _m_ is the group size of observations removed from the sample for jackknifing (Busing et al. 1999). See our corresponding manuscript "_Genomic signals of admixture as well as reinforcement between two sister species of European sepsid flies in sympatry vs. parapatry_" (Giessen _et al._ 2020) for more details on the approach and differences to other methods

## Workflow

The ususal starting point for analyses of geneflow with ABBABABA-4AF is the generation of an Allele frequency matrix file. This can be done, for example, by converting a SNP file in the SYNC file format (Kofler _et al._ 2011) to major allele frequencies using the provided python script "SYNC2AF.py". Files in the SYNC file format are tab-delimited, where the first three columns are "Chromosome", "Position" and "Reference Allele", respectively. These are followed by columns of allele counts for one or multiple synchronized population samples. The allele counts are encoded in the form: "A:T:C:G:N:Del". For example, a population with 50% A/G polymorphism at a read depth 20, would be encoded as: "10:0:0:10:0:0".

A typical command line for SYNC2AF.py looks like this:

```bash
python3 SYNC2AF.py --sync input.sync.gz | gzip > output.af.gz
```

The output file is tab-delimited where the first three columns are "Chromosome", "Position" and "Major/Minor alleles", respectively. These are followed by the allele frequencies of the major allele for a given position for each of the populations in the sync file. Note, that only the two most common alleles are considered and counts of other alleles at tri- or tetra-allelic position are ignored.

This file can then be used as the input for the ABBA-BABA test. ABBABABA-4AF.py requires the following input files and parameters:

-   AlleleFrequencies: The 2D matrix of allele frequencies as described above
-   SNPs: The number of SNPs in each window per contig or chromosome
-   Order: The column position of the four populations (H1,H2,H3 and H4) considered for the ABBA-BABA test. Note, that Python uses zero-based indexing. Thus, if you want to use the populations in the 4<sup>th</sup>, 8<sup>th</sup>, 9<sup>th</sup> and 10<sup>th</sup> columns as H1, H2, H3 and H4, you will need to provide a comma-separated list 3,7,8,9 as input.
-   Output: You need to provide an output-prefix which will then be used to generate two output files with the endings "\_GenomeWide.D" and "\_XXXSNPs.D", where XXX is replaced by the number of SNPs provided above.

A typical command line for ABBABABA-4AF.py (with focal populations in the 4<sup>th</sup>, 8<sup>th</sup>, 9<sup>th</sup> and 10<sup>th</sup> columns of the input) looks like this:

```bash
python3 ABBABABA-4AF.py --AF output.af.gz --SNPs 500 --Order 3,7,8,9 --Output output-file
```

### Description of Output Files

#### output-file_500SNPs.D

This file is the output containing window-wise estimates of _D_<sub>S</sub> and _D_<sub>D</sub>. Here the columns are:

-   Chromosome
-   Start Position
-   Average Position
-   End Position
-   Window Size in basepairs
-   Number of SNPs in Window
-   Soraggi's _D_
-   Durand's _D_

#### output-file_genome-wide.D

This file is the output containing genome-wide estimates of _D_<sub>S</sub> and _D_<sub>D</sub> as well as jackknifing variance and _z_-scores. Here the columns are:

-   Analysis Method (Soraggi or Durand)
-   Genome-wide estimate of _D_
-   Jackknife estimate of _D_
-   Jackknife estimate of the variance of _D_
-   _Z_-score

## References

Busing, F.M.T.A., Meijer, E. & Leeden, R.V.D. 1999. Delete-m Jackknife for Unequal m. Statistics and Computing 9: 3–8.

Durand, E.Y., Patterson, N., Reich, D. & Slatkin, M. 2011. Testing for Ancient Admixture between Closely Related Populations. Mol Biol Evol 28: 2239–2252.

Green, R.E., Krause, J., Briggs, A.W., Maricic, T., Stenzel, U., Kircher, M., et al. 2010. A Draft Sequence of the Neandertal Genome. Science 328: 710–722.

Kofler, R., Pandey, R.V. & Schlotterer, C. 2011. PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq). Bioinformatics 27: 3435–3436.

Soraggi, S., Wiuf, C. & Albrechtsen, A. 2018. Powerful Inference with the D-Statistic on Low-Coverage Whole-Genome Data. G3: Genes, Genomes, Genetics 8: 551–566.
