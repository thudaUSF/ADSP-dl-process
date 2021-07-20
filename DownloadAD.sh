#!/bin/zsh
#SBATCH --job-name=DownloadAD
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8GB

#module purge
#module add apps/samtools/1.3.1

#BESURE TO CHANGE i (line 14), end (line 15), and ESPECIALLY URL (line 36)
input="/home/ec2-user/manifest.txt" #
DownloadFolder="/home/ec2-user/bams/crams" #Where downloads go
Craifolder="/home/ec2-user/crai/"
Workdir="/home/ec2-user/"

i=1 #CHANGE THIS; Start from 4 (including)
end=$i+10 #CHANGE THIS; Count from 4 - 12 (Does 10 files)
n=0
echo "Deleted:"
find ${DownloadFolder} -maxdepth 1 -name "*.bam" -type 'f' -size -100k -print -delete #delete empty files only in the current directory
normalsize=100000

cd $Craifolder #samtools afaik can only read the crai if it's in the location (see URL)

while read -r line; do
	
	((n++)) #Counter, stops when it does amount of files 
	if (($n == $end)); then 
		break
	fi

	if (($n>=$i)); then #If the filename is within the selected choices; Double parantheses are used to use "normal expressions"
		IDcram=$(cut -d',' -f1 <<< "$line")
		snd=$(cut -d',' -f2 <<< "$line")
		ID="${IDcram}_vcpa1.1"
		echo "Downloading/Checking: ${ID} - ${n}"
		#ID=${ID//$'\r'} get rid of carriage return, which resets things after a variable; fixed bc manifest was written in DOS
		
		#samtools merge /home/ec2-user/bams/crams/mergedfolder/"${ID}.bam" /home/ec2-user/bams/crams/${ID}*.bam #for merging files into one
		
		#if grep -Fxq "${ID}" "${Workdir}/completed.txt"; then #if every one of them is completed, don't bother doing anything else
		#	echo "Above file was already completed, continuing" #should be continue 
		#	continue
		#fi
		
		#Download crai file, index file needed for any of them
		if [ ! -f "${Craifolder}/${ID}.cram.crai" ]; then 
			string="s3://wanglab-dss-share/distribution/adsp/cram/${snd}/${ID}.cram.crai"
			echo $string
			aws s3 cp --request-payer requester $string "${Craifolder}/${ID}.cram.crai"
		fi
		
		URLhts=`python3 /home/ec2-user/Scripts/urlget.py $snd $ID` #result ends up like - ./htsfile -h 'https://wanglab-dss-share.s3.amazonaws.com/distribution/adsp/cram/filename.cram?AWSAccessKeyId=ID&Signature=Signature&x-amz-request-payer=requester&Expires=1605118813'
		URL=${URLhts:14} #cut off "./htsfile -h '"
		URL=${URL::-1} #cut off " end '"
		URL="${URL}&local=/${ID}.cram"
		
		if [ ! -f "/home/ec2-user/bams/HLA/${ID}HLA.bam" ]; then 
			samtools view --reference /home/ec2-user/bams/Homo_sapiens_assembly38.fasta $URL chr6:29844528-33100696 -b > "/home/ec2-user/bams/HLA/${ID}HLA.bam"
			echo "HLA complete"
			filesizeHLA=$(stat -c%s "/home/ec2-user/bams/HLA/${ID}HLA.bam")
		fi
		if [ ! -f "${DownloadFolder}/${ID}TRA.bam" ]; then 
			samtools view --reference /home/ec2-user/bams/Homo_sapiens_assembly38.fasta $URL chr14:21621904-22552132  -b > "${DownloadFolder}/${ID}TRA.bam"
			samtools index "${DownloadFolder}/${ID}TRA.bam"
			echo "TRA complete"
			filesizeTRA=$(stat -c%s "${DownloadFolder}/${ID}TRA.bam")
		fi
		if [ ! -f "${DownloadFolder}/${ID}TRD.bam" ]; then 
			samtools view -o "${DownloadFolder}/${ID}TRD.bam" -b "${DownloadFolder}/${ID}TRA.bam" chr14:22422546-22466577 #get TRD from within TRA
			echo "TRD complete"
			filesizeTRD=$(stat -c%s "${DownloadFolder}/${ID}TRD.bam")
		fi
		if [ ! -f "${DownloadFolder}/${ID}TRB.bam" ]; then 
			samtools view --reference /home/ec2-user/bams/Homo_sapiens_assembly38.fasta $URL chr7:142299011-142813287 -b > "${DownloadFolder}/${ID}TRB.bam"
			echo "TRB complete"
			filesizeTRB=$(stat -c%s "${DownloadFolder}/${ID}TRB.bam")
		fi
		if [ ! -f "${DownloadFolder}/${ID}TRG.bam" ]; then 
			samtools view --reference /home/ec2-user/bams/Homo_sapiens_assembly38.fasta $URL chr7:38240024-38368055 -b > "${DownloadFolder}/${ID}TRG.bam"
			echo "TRG complete"
			filesizeTRG=$(stat -c%s "${DownloadFolder}/${ID}TRG.bam")
		fi
		if [ ! -f "${DownloadFolder}/${ID}IGH.bam" ]; then 
			samtools view --reference /home/ec2-user/bams/Homo_sapiens_assembly38.fasta $URL chr14:105586437-106879844  -b > "${DownloadFolder}/${ID}IGH.bam"
			echo "IGH complete"
			filesizeIGH=$(stat -c%s "${DownloadFolder}/${ID}IGH.bam")
		fi
		if [ ! -f "${DownloadFolder}/${ID}IGK.bam" ]; then 
			samtools view --reference /home/ec2-user/bams/Homo_sapiens_assembly38.fasta $URL chr2:88857361-90235368 -b > "${DownloadFolder}/${ID}IGK.bam" 
			echo "IGK complete"
			filesizeIGK=$(stat -c%s "${DownloadFolder}/${ID}IGK.bam")
		fi
		if [ ! -f "${DownloadFolder}/${ID}IGL.bam" ]; then 
			samtools view --reference /home/ec2-user/bams/Homo_sapiens_assembly38.fasta $URL chr22:22026076-22922913 -b > "${DownloadFolder}/${ID}IGL.bam"
			echo "IGL complete"
			filesizeIGL=$(stat -c%s "${DownloadFolder}/${ID}IGL.bam")
		fi

		
		#Unmapped region		
		#if [ ! -f "${DownloadFolder}/${ID}unmapped.bam" ]; then 
			#samtools view --reference /home/ec2-user/bams/Homo_sapiens_assembly38.fasta $URL '*' -b > "${DownloadFolder}/${ID}unmapped.bam"
			#samtools view --reference /home/ec2-user/bams/Homo_sapiens_assembly38.fasta -b -f 4 $URL > "${DownloadFolder}/${ID}unmapped.bam"
			#echo "Unmapped complete"
		#fi 
		
		#echo "${filesizeHLA}, ${filesizeTRA}, ${filesizeTRB}, ${filesizeTRG}, ${filesizeIGH}, ${filesizeIGK}, ${filesizeIGL}"
		if (( filesizeHLA > normalsize)) && (( filesizeTRA > normalsize)) && (( filesizeTRD > normalsize)) && (( filesizeTRB > normalsize)) && (( filesizeTRG > normalsize)) && (( filesizeIGH > normalsize)) && (( filesizeIGK > normalsize)) && (( filesizeIGL > normalsize)); then 
			echo -e "${ID}" >> "${Workdir}/completed.txt"
		fi
		
	fi

done < $input