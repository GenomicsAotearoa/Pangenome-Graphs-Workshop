#!/usr/bin/bash
#Handling command line options
SHORT=r:,s:,i:,d:,l:,o:h
LONG=ref:,snp:,indel:,depth:,length:,output:,help
OPTS=$(getopt -a -n sim_vc_compare  --options $SHORT --longoptions $LONG -- "$@")
eval set -- "$OPTS"
while :
do
  case "$1" in
    -r | --ref )
      ref="$2"
      shift 2
      ;;
    -s | --snp )
      snp="$2"
      shift 2
      ;;
    -i | --indel )
      indel="$2"
      shift 2
      ;;
    -d | --depth )
      depth="$2"
      shift 2
      ;;
    -l | --length )
      length="$2"
      shift 2
      ;;
    -o | --output )
      output="$2"
      shift 2
      ;;
    -h | --help)
      echo "Program : sim_vc_compare"
      echo "Version : 1.0"
      echo "Contact : fathima.nuzla.ismail@gmail.com"
      echo "Usage   : sim_vc_compare.sh [options]"
      echo "Options :"
      echo "-r | --ref STR reference sequence file"
      echo "-s | --snp INT Number of SNPs to simulate"
      echo "-i | --indel INT Number of INDELs to simulate"
      echo "-d | --depth INT Cover depth of the reads"
      echo "-l | --length INT of a read"
      echo "-o | --output STR Output folder name"
      echo "-h | --help Display this help message" 
      exit 2
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      ;;
  esac
done

# exit when any command fails
set -e

if command -v module &> /dev/null
then
echo ">>> Loading required modules ...";
#Load required modules with specific versions
   module purge
   module load BCFtools/1.9-GCC-7.4.0
   module load SAMtools/1.9-GCC-7.4.0
   module load BWA/0.7.17-GCC-9.2.0
   module load wgsim/20111017-GCC-11.3.0
fi;

echo ">>> Checking for reference file..."
if [ -z "$ref" ] || [ ! -f "$ref" ]; then
    echo "Reference file ${ref} does not exists or not specifies by --ref <filename> !"
    exit 1;
fi

#Check arguments --snp and --indel
if [ -z $snp ] && [ -z $indel ]; then
   echo "One of --snp or --indel must be specified !"
   exit 1;
fi

echo "Checking for simuG.pl file..."
if [ -f "simuG.pl" ]; then
    echo "simuG.pl script is available !"
else 
    echo "simuG.pl script is inot available !"
    echo ">>> Downloading it from git ..."
    wget https://raw.githubusercontent.com/yjx1217/simuG/master/simuG.pl;
fi

sleep 2;

echo ">>> Finding the length of the reference file ..."
samtools faidx $ref;
reflen=$(head -1 ${ref}.fai | cut -f 2)
echo "Length is ${reflen}"

echo ">>> Verifying snp and indel counts ..."
if [ -z $snp ]; then
   snp=0
fi

if [ -z $indel ]; then
   indel=0
fi

re='^[0-9]+$'
if ! [[ $snp =~ $re ]] ; then
   echo "--snp must be an integer !"; exit 1;
fi

if ! [[ $indel =~ $re ]] ; then
   echo "--indel must be an integer !"; exit 1;
fi

if [ "$snp" -gt $(($reflen/2)) ]; then
  echo "--snp must be less than half of the reference sequence length !";
fi

if [ "$indel" -gt $(($reflen/2)) ]; then
  echo "--indel must be less than half of the reference sequence length !";
fi

echo ">>> Checking for the output directory ...";
if [ -z $output ]; then
    output="output"
    echo "No output directory is specified. using default directory \"output\"."
fi
if [ -d $output ]; then
    echo "The directory \"${output}\" exists. Files will be overriden !";
else 
    echo "Creating the directory \"${output}\" ..."
    mkdir $output;
fi

echo ">>> Changing the directory to $(pwd) ..."
cd $output;
echo ">>> Simulating the new sequence using simuG ...";
simprefix="sim_snp_${snp}_indel_${indel}";
perl ../simuG.pl -refseq "../${ref}" -snp_count $snp -indel_count $indel -prefix $simprefix;

sleep 5;

echo ">>> Renaming chromosome name of the simulated sequence ..."
sed -i '/>/ s/\(.*\)/\1_'${simprefix}'/' "${simprefix}.simseq.genome.fa"

echo ">>> bgzip and merging output vcf files ..."
bgzip -f "${simprefix}.refseq2simseq.SNP.vcf"
bcftools index "${simprefix}.refseq2simseq.SNP.vcf.gz"
bgzip -f "${simprefix}.refseq2simseq.INDEL.vcf"
bcftools index "${simprefix}.refseq2simseq.INDEL.vcf.gz"

sleep 5;

bcftools merge "${simprefix}.refseq2simseq.SNP.vcf.gz" "${simprefix}.refseq2simseq.INDEL.vcf.gz" -O z -o "${simprefix}.gt.vcf.gz"
bcftools index "${simprefix}.gt.vcf.gz"

echo ">>> Calculating number of simulated reads required ..."
if [ -z $depth ]; then
   depth=30;
   echo "No coverage depth is specified. Assigning default 30"; 
fi
if ! [[ $depth =~ $re ]] ; then
   echo "--depth must be an integer !"; exit 1;
fi

if [ -z $length ]; then
   length=100;
   echo "No reads length is specified. Assigning default 100";
fi
if ! [[ $length =~ $re ]] ; then
   echo "--length must be an integer !"; exit 1;
fi
echo ">>> Calculating number of reads required ..."
reads=$(($reflen*$depth/$length))
echo "Need ${reads} reads for ${depth} coverage depth with ${length} reads";
echo ">>> Simulating reads ..."
wgsim -N ${reads} -1 ${length} -2 ${length} "${simprefix}.simseq.genome.fa" "${simprefix}.read1.fq" "${simprefix}.read2.fq";

sleep 5;

echo ">>> Aligning reads to the reference \"${ref}\" and generating the VCF file using bwa mem ...";
bwa index "../${ref}"
bwa mem -R "@RG\tID:${simprefix}\tSM:${simprefix}\tLB:L1" "../${ref}" "${simprefix}.read1.fq" "${simprefix}.read2.fq" > "${simprefix}.sam"
samtools view -bS "${simprefix}.sam" | samtools sort - > "${simprefix}.bam"
bcftools mpileup -Ou -f "../${ref}" "${simprefix}.bam" | bcftools call -vmO z -o "${simprefix}.bwa.vcf.gz"
bcftools index "${simprefix}.bwa.vcf.gz"

sleep 5;

echo ">>> Comparing the ground truth vcf and simulated vcf using bcftools isec ..."
bcftools isec -c none -p compare "${simprefix}.gt.vcf.gz" "${simprefix}.bwa.vcf.gz"

sleep 5;

echo "";
bold=$(tput bold)
normal=$(tput sgr0)
echo "${bold}>>> Generating stats ...${normal}";
sim_snp=$(bcftools stats "${simprefix}.bwa.vcf.gz" | grep "number of SNPs:" | cut -f 4)
sim_indel=$(bcftools stats "${simprefix}.bwa.vcf.gz" | grep "number of indels:" | cut -f 4)
sim_u_snp=$(bcftools stats compare/0001.vcf | grep "number of SNPs:" | cut -f 4)
sim_u_indel=$(bcftools stats compare/0001.vcf | grep "number of indels:" | cut -f 4)
sim_c_snp=$(bcftools stats compare/0002.vcf | grep "number of SNPs:" | cut -f 4)
sim_c_indel=$(bcftools stats compare/0002.vcf | grep "number of indels:" | cut -f 4)

sleep 1;

tp=$(($sim_c_snp+$sim_c_indel))
fp=$(($sim_u_snp+$sim_u_indel))
tn=$(($reflen-$snp-$indel-$sim_u_snp-$sim_u_indel))
fn=$(($snp+$indel-$tp))
sensityvity=$(bc <<< "scale=4; (${tp}*100/(${tp}+${fn}));")
specificity=$(bc <<< "scale=4; (${tn}*100/(${tn}+${fp}));")
echo   "--------------------------------------------------"
echo   "|  ${bold}REPORT${normal}                                        |"
echo   "--------------------------------------------------"
printf "|  Ground Truth SNPs                = %'10d |\n" ${snp}                  
printf "|  Ground Truth INDELs              = %'10d |\n" ${indel}                
printf "|  Identified SNPs in Simulation    = %'10d |\n" ${sim_snp}              
printf "|  Identified INDELs in Simulation  = %'10d |\n" ${sim_indel}            
printf "|  SNPs Private to Simulation       = %'10d |\n" ${sim_u_snp}            
printf "|  INDELs Private to Simulation     = %'10d |\n" ${sim_u_indel}          
printf "|  Exact Matched SNPs               = %'10d |\n" ${sim_c_snp}            
printf "|  Exact Matched INDELs             = %'10d |\n" ${sim_c_indel}            
printf "|  True Positive (TP)               = %'10d |\n" ${tp}
printf "|  False Positive (FP)              = %'10d |\n" ${fp}
printf "|  True Negative (TN)               = %'10d |\n" ${tn}
printf "|  False Negative (FN)              = %'10d |\n" ${fn}
echo   "|------------------------------------------------"
printf "|  ${bold}Sensitivity                      = %'9.4f%%${normal} |\n" ${sensityvity} 
printf "|  ${bold}Specificity                      = %'9.4f%%${normal} |\n" ${specificity} 
echo "--------------------------------------------------"
echo "End of the program !"
 
