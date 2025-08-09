[genome_stats.md](https://github.com/user-attachments/files/21678285/genome_stats.md)[genome_stats.md](https://github.com/user-attachments/files/21678283/genome_stats.md)
# BioinformaticsProject
Our project for CSS383

<ins>Project 1 Proposal: Are Genes in Symbiotic Bacteria Smaller Than in Free-Living Bacteria?<ins>

Symbiotic bacteria, especially those that live inside host cells, often go through what’s called genome reduction. Basically, they lose a lot of the DNA they don’t need anymore because their environment is stable and predictable. Scientists have known for a while that these bacteria tend to have much smaller genomes overall, but one question that hasn’t been explored as much is whether their individual genes are also getting smaller.
This project focuses on comparing gene lengths between symbiotic bacteria and closely related free-living ones to see if the symbiotic lifestyle leads to shorter genes. If it does, it might tell us more about how living inside a host affects not just what genes bacteria have, but how those genes are structured.

<ins>Background studies:<ins>

McCutcheon & Moran (2011) – Found that highly reduced symbiont genomes (e.g., Buchnera, Carsonella) pack genes closely, sometimes with overlaps. While gene length wasn’t directly measured, the patterns suggest shorter genes may be part of the process.

Giovannoni et al. (2005) – Studied free-living marine bacteria (Prochlorococcus, Pelagibacter ubique) with small genomes caused by nutrient limitation and strong selection for efficiency, not by relaxed selection.

Moran & Bennett (2014) – Review of ultra-reduced insect symbiont genomes; synthesizes mechanisms of genome erosion (gene loss, drift, deletion bias) and highlights record-small genomes.

Wang et al. (2024) – Compared symbiotic Fluviibacter to free-living relatives. Found symbiont genes averaged 545 bp vs. 943 bp in free-living counterparts. Also observed higher dN/dS ratios in symbionts, indicating genetic drift.

Boscaro et al. (2013) – Examined Polynucleobacter necessarius, showing major gene loss in symbionts, suggesting simplification.

Moran & Bennett (2014) – Reviewed symbiont genome reduction, noting gene loss, drift, and shorter/overlapping genes.

<ins>Contrast:<ins>

Symbionts may have shorter genes due to relaxed selective pressures and deletion bias.

Free-living bacteria may have shorter genes due to selection for metabolic efficiency.

<ins>References:<ins>

Boscaro, V. et al. (2013). Polynucleobacter necessarius, a model for genome reduction in both free-living and symbiotic bacteria. PNAS, 110(46), 18590-18595. https://www.pnas.org/doi/pdf/10.1073/pnas.1316687110

McCutcheon, J.P., & Moran, N.A. (2011). Extreme genome reduction in symbiotic bacteria. Nature Reviews Microbiology, 10(1), 13-26. https://www.nature.com/articles/nrmicro2670

Giovannoni, S.J. et al. (2005). Genome streamlining in a cosmopolitan oceanic bacterium. Science, 309(5738), 1242-1245. https://www.science.org/doi/10.1126/science.1114057

Bennett, G.M., & Moran, N.A. (2013). Small, smaller, smallest: the origins and evolution of ancient dual symbioses in a phloem-feeding insect. Genome Biol Evol, 5(9), 1675-1688. https://academic.oup.com/gbe/article/5/9/1675/555845

Wang, R. et al. (2024). Comparative genomic analysis of symbiotic and free-living Fluviibacter phosphoraccumulans strains. Appl Environ Microbiol, 90(3), e01900-23. https://journals.asm.org/doi/epub/10.1128/aem.01900-23

Moran, N.A., & Wernegreen, J.J. (2000). Lifestyle evolution in symbiotic bacteria: insights from genomics. Trends Ecol Evol, 15(8), 321-326. https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=2bcff42d8c43855bf191da21e817ef8a05c756b1

Methanosaeta Concilli Genome - https://www.ncbi.nlm.nih.gov/nuccore/CP002565.1?report=genbank&log$=seqview

Buchnera aphidicola str. APS (Acyrthosiphon pisum) APS DNA - https://www.ncbi.nlm.nih.gov/nuccore/AP036055.1?report=genbank&log$=seqview

Salmonella enterica strain no75 chromosome, complete genome - https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP075372.1

Escherichia coli str. K-12 substr. MG1655 - https://www.ncbi.nlm.nih.gov/nuccore/U00096.3

Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis DNA, complete genome - https://www.ncbi.nlm.nih.gov/nuccore/BA000021.3?report=genbank&log$=seqview

<ins>Bacteria Being Compared<ins>

[genome_stats.md](https://github.com/user-attachments/files/21696129/genome_stats.md)
| Bacterium | FASTA Path | Genome Size (bp) | GC Content (%) |
|---|---|---:|---:|
| Buchnera aphidicola (symb) | `Buchnera aphidicola str. APS (Acyrthosiphon pisum) APS DNA.fasta` | 641802 | 26.30 |
| E. coli K-12 MG1655 (free) | `Escherichia coli str. K-12 substr. MG1655.fasta` | 4641652 | 50.79 |
| Wigglesworthia glossinidia (symb) | `wigglesworthia_g_sequence.fasta` | 697724 | 22.48 |
| Salmonella enterica no75 (symb) | `salmonella_sequence.fasta` | 5032348 | 52.16 |
| Chlamydophila pneumoniae TW-183 (symb) | `c_pneumoniae.fasta` | 1225935 | 40.58 |
| Isosphaera pallida ATCC 43644 (free) | `i_pallida.fasta` | 5472964 | 62.44 |


<ins>Comparison Pairs:<ins>

Symbiotic: Buchnera aphidicola str. APS
Free-living counterpart: Escherichia coli str. K-12 MG1655

Symbiotic: Wigglesworthia glossinidia
Free-living counterpart: Salmonella enterica strain no75

Symbiotic: Salmonella enterica strain no75
Free-living counterpart: Methanosaeta concilii GP-6

<ins>Hypothesis:<ins>

Symbiotic bacteria have shorter genes than their free-living counterparts due to relaxed selective pressures and deletion bias in stable host environments.

<ins>Predictions:<ins>

Symbiotic genes will be significantly shorter on average.

The difference will be most pronounced in genes related to regulation and signaling.

Housekeeping genes will show smaller differences.

<ins>Approach:<ins>

Extract gene lengths from GFF files using Python (Biopython, Pandas).

Classify genes by bacterial lifestyle (symbiotic or free-living).

Statistical analysis — Wilcoxon signed-rank test for paired distributions.

Visualization — Matplotlib plots of gene length distributions and possible genome size trends.

<ins>Group Members:<ins>

Angela Monsegue, Destiny Rebman, Arham Saed, Katya Sizov
