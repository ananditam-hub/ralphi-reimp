# ralphi-reimp
This repository attempts to reimplement and reproduce the results obtained in the paper titled: 'ralphi: a deep reinforcement learning framework for haplotype assembly'. <br/><br/>
We implemented RALPHI on the NA12878 human genome using Illumina sequencing data at 5× and 30× coverage. Use this link to download the .vcf and .cram files corresponding to NA12878: [International Genome Sample Resource (NA12878.final.cram)](https://www.internationalgenome.org/data-portal/sample/NA12878). <br/>
Unlike WhatsHap, Ralphi isn't compatible with .cram files. The following commands of samtools can be used to convert the .cram file to .bam file and index the same: `samtools view -b -o output.bam input.cram` and `samtools index input.bam`. The latter results in a .bam.bai file which is the index file and should be placed in the same directory as the .bam. <br/>
The above procedure can be time and memory heavy since the file is **10 GB** huge, thus it is advised to use suitable resources when using this phasing tool. <br/>
For execution, update the illumina/ont.yaml files in the n_config directory, based on the source of reads, with suitable path labels to .bam, .vcf, .fasta (optional for Illumina reads) files. Thereafter you can either use the terminal and execute: `conda activate ralphi-env` and `python3 phase.py --config \path\to\yaml` or the shell script:`shell_script/ralphi_exe.sh`. You can also follow the instructions on the github repository: [https://github.com/PopicLab/ralphi](https://github.com/PopicLab/ralphi)
<br/>
## Result Summary
On implementing ralphi on chromosomes 4, 14 and 22 (chosen at random), of NA12878 with 5x coverage, without reference data, we get the following plot:<br/>
<img src="results/switch_rate_5x.png" alt="5x Coverage" width="350"/><br/>
And for chromosomes 1, 14 and 22 with 30x coverage (we leave out LongPhase here due to longer execution time, will update soon): <br/>
<img src="results/switch_rate_30x.png" alt="30x Coverage" width="350"/><br/>
Phasing quality, as measured by switch rate, improves consistently from low to high coverage across all tools. Ralphi makes significant improvements in achieving accuracy from its previous formulation as in RefHap based on maximum fragment cut formulation. Its results can be seen comparable to that of HapCUT2, a tool based on Minimum Error Correction (MEC) formulation, another state-of-the-art technique. <br/><br/>
<u>**Time Comparison**</u> <br/> For chr22 of 30x coverage NA12878 human genome the following is the time taken by each of the four* tools, without using a reference genome: <br/>
**WhatsHap:** 3m47.923s <br/>
**Ralphi:** 3m17.995s <br/>
**HapCUT2:** 2m27.970s (for all chromosomes).



## Other Tools Used:
Samtools: [https://github.com/samtools/samtools.git](https://github.com/samtools/samtools.git) <br/>
Bcftools: [https://github.com/samtools/bcftools.git](https://github.com/samtools/bcftools.git) <br/>
WhatsHap: [https://github.com/whatshap/whatshap/tree/bb7ccfffc655072451d642b4eea9661f96b345af](https://github.com/whatshap/whatshap/tree/bb7ccfffc655072451d642b4eea9661f96b345af) <br/>
Longphase:[https://github.com/twolinin/longphase](https://github.com/twolinin/longphase), or you can use the conda package, which is easier to use: run `conda install bioconda::longphase` in terminal.<br/>
HapCUT2: [https://github.com/vibansal/HapCUT2](https://github.com/vibansal/HapCUT2), or you can use the conda package, which is easier to use: run `conda install bioconda::hapcut2` in terminal.<br/>
## Citation
`Battistella, E., Maheshwari, A., Ekim, B., Berger, B., & Popic, V. (2025). ralphi: a deep reinforcement learning framework for haplotype assembly. bioRxiv : the preprint server for biology, 2025.02.17.638151. https://doi.org/10.1101/2025.02.17.638151`
