# ralphi-reimp
This repository attempts to reimplement and reproduce the results obtained in the paper titled: 'ralphi: a deep reinforcement learning framework for haplotype assembly'. <br/><br/>
Use this link to download the .cram file corresponding to NA12878: ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram <br/>
Unlike Whatshap, Ralphi isn't compatible with .cram files. The following commands of samtools can be used to convert the .cram file to .bam file and index the same: `samtools view -b -o output.bam input.cram` and `samtools index input.bam`. The latter results in a .bam.bai file which is the index file and should be placed in the same directory as the .bam. <br/>
