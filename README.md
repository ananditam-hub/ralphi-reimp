# ralphi-reimp
This repository attempts to reimplement and reproduce the results obtained in the paper titled: 'ralphi: a deep reinforcement learning framework for haplotype assembly'. <br/><br/>
Use this link to download the .cram file corresponding to NA12878: [International Genome Sample Resource](ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram). <br/>
Unlike WhatsHap, Ralphi isn't compatible with .cram files. The following commands of samtools can be used to convert the .cram file to .bam file and index the same: `samtools view -b -o output.bam input.cram` and `samtools index input.bam`. The latter results in a .bam.bai file which is the index file and should be placed in the same directory as the .bam. <br/>
The above procedure can be time and memory heavy since the file is **10 GB** huge, thus it is advised to used suitable resources when using this phasing tool.
<br/>
## Other Tools Used:
Samtools: [https://github.com/samtools/samtools.git](https://github.com/samtools/samtools.git) <br/>
Bcftools: [https://github.com/samtools/bcftools.git](https://github.com/samtools/bcftools.git) <br/>
WhatsHap: [https://github.com/whatshap/whatshap/tree/bb7ccfffc655072451d642b4eea9661f96b345af](https://github.com/whatshap/whatshap/tree/bb7ccfffc655072451d642b4eea9661f96b345af) <br/>
Longphase:[https://github.com/twolinin/longphase](https://github.com/twolinin/longphase), or you can use the conda package, which is easier to use: run `conda install bioconda::longphase` in terminal.<br/>
HapCUT2: [https://github.com/vibansal/HapCUT2](https://github.com/vibansal/HapCUT2), or you can use the conda package, which is easier to use: run `conda install bioconda::hapcut2` in terminal.<br/>
## Citation
Battistella, E., Maheshwari, A., Ekim, B., Berger, B., & Popic, V. (2025). ralphi: a deep reinforcement learning framework for haplotype assembly. bioRxiv : the preprint server for biology, 2025.02.17.638151. https://doi.org/10.1101/2025.02.17.638151
