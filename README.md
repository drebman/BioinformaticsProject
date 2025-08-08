[genome_stats.md](https://github.com/user-attachments/files/21678285/genome_stats.md)[genome_stats.md](https://github.com/user-attachments/files/21678283/genome_stats.md)# BioinformaticsProject
Our project for CSS383

Project 1 Proposal: Are Genes in Symbiotic Bacteria Smaller Than in Free-Living Bacteria?

Symbiotic bacteria, especially those that live inside host cells, often go through what’s called genome reduction. Basically, they lose a lot of the DNA they don’t need anymore because their environment is stable and predictable. Scientists have known for a while that these bacteria tend to have much smaller genomes overall, but one question that hasn’t been explored as much is whether their individual genes are also getting smaller.
This project focuses on comparing gene lengths between symbiotic bacteria and closely related free-living ones to see if the symbiotic lifestyle leads to shorter genes. If it does, it might tell us more about how living inside a host affects not just what genes bacteria have, but how those genes are structured.

Bacteria Being Compared

| Bacterium | FASTA Path | Genome Size (bp) | GC Content (%) |
|---|---|---:|---:|
| Buchnera aphidicola (symb) | `Buchnera aphidicola str. APS (Acyrthosiphon pisum) APS DNA.fasta` | 641802 | 26.30 |
| Wigglesworthia glossinidia (symb) | `wigglesworthia_g_sequence.fasta` | 697724 | 22.48 |
| Salmonella enterica no75 (symb) | `salmonella_sequence.fasta` | 5032348 | 52.16 |
| E. coli K-12 MG1655 (free) | `Escherichia coli str. K-12 substr. MG1655.fasta` | 4641652 | 50.79 |
| Methanosaeta concilii GP-6 (free) | `Methanosaeta concilii GP-6.fasta` | 3008626 | 51.03 |

Comparison Pairs
Symbiotic: Buchnera aphidicola str. APS
Free-living counterpart: Escherichia coli str. K-12 MG1655

Symbiotic: Wigglesworthia glossinidia
Free-living counterpart: Salmonella enterica strain no75

Symbiotic: Salmonella enterica strain no75
Free-living counterpart: Methanosaeta concilii GP-6
