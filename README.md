# ADSP-dl-process
Scripts to download exome slices from ADSP (Alzheimer's Disease Sequencing Project) (Requires bucket connection to ADSP site) + processing for t-test to compare physchem 
# For downloading:
- Need to install prerequisites, samtools, go, fuse-utils, goofys, then enter:
  - export GOPATH=$HOME/work
  - $GOPATH/bin/goofys adsp-exomes /home/ec2-user/bams 
  - ^bucket named adsp-exomes goes to the bams folder, allows for large dl storage
- Download reference hg38 file from broad-institute 
- Connect to Filezilla using ec2 token, URL from EC2 instance
- If htslib has issues:
  - export LD_LIBRARY_PATH=/home/ec2-user/htslib-1.11
- ADdownload.sh depends on urlget.py, run ADdownload.sh
