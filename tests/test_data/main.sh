# simulate BAM and FASTQ data from reference fasta and vcf file

python2 /Users/charlesmurphy/Desktop/tools/neat-genreads/genReads.py \
  -r ref.fasta \
  -R 50 \
  -o test_sample1 \
  -c 30 \
  -p 12 \
  -v test.vcf \
  --bam \
  -M 0 \
  -E 0 \
  --pe 120 5
~/Desktop/tools/samtools-1.10/samtools index test_sample1_golden.bam

python2 /Users/charlesmurphy/Desktop/tools/neat-genreads/genReads.py \
  -r ref.fasta \
  -R 50 \
  -o test_sample2 \
  -c 30 \
  -p 12 \
  -v test.vcf \
  --bam \
  -M 0 \
  -E 0 \
  --pe 120 5
~/Desktop/tools/samtools-1.10/samtools index test_sample2_golden.bam

# python biometrics/cli.py extract -sb ./tests/test_data/test_golden.bam -st tumor -ss male -sp P1 -sn test --vcf tests/test_data/test.vcf -db . --fafile ./tests/test_data/ref.fasta --overwrite
#
# python biometrics/cli.py minor -sb ./tests/test_data/test_golden.bam -st tumor -ss male -sp P1 -sn test --vcf tests/test_data/test.vcf -db . --fafile ./tests/test_data/ref.fasta
#
# python biometrics/cli.py genotype \
#   -sb ./tests/test_data/test_sample1_golden.bam ./tests/test_data/test_sample2_golden.bam \
#   -st tumor tumor -ss male male -sg patien1 patient1 -sn test_sample1 test_sample2 \
#   --vcf tests/test_data/test.vcf -db . --fafile ./tests/test_data/ref.fasta
#
#
# python biometrics/cli.py extract \
#   -sb /Users/charlesmurphy/Desktop/mskcc-analyses/200608_compare_qc_tools/manually_count_bases/C-52YNHF-N001-d_cl_aln_srt_MD_IR_FX_BR.bam \
#   -st tumor -ss male -sg C-52YNHF -sn C-52YNHF-N001-d \
#   --vcf /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/fingerprinting_snps.vcf \
#   -db /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/ \
#   --fafile /Users/charlesmurphy/Desktop/data/ref/hg19/Homo_sapiens_assembly19.fasta \
#   --overwrite
#
biometrics extract \
  -sb /Users/charlesmurphy/Desktop/mskcc-analyses/200608_compare_qc_tools/manually_count_bases/C-52YNHF-N001-d_cl_aln_srt_MD_IR_FX_BR.bam \
  -st tumor -ss male -sg C-52YNHF -sn C-52YNHF-N002-d \
  --vcf /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/fingerprinting_snps.vcf \
  --bed /Users/charlesmurphy/Desktop/data/innovation/resources/MSK-ACCESS-v1.0/MSK-ACCESS-v1_0-probe-A.sorted.bed \
  -db /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/db/ \
  --fafile /Users/charlesmurphy/Desktop/data/ref/hg19/Homo_sapiens_assembly19.fasta \
  --overwrite
#
biometrics genotype \
  -sb \
  /Users/charlesmurphy/Desktop/mskcc-analyses/200608_compare_qc_tools/manually_count_bases/C-52YNHF-N001-d_cl_aln_srt_MD_IR_FX_BR.bam \
  /Users/charlesmurphy/Desktop/mskcc-analyses/200608_compare_qc_tools/manually_count_bases/C-52YNHF-N001-d_cl_aln_srt_MD_IR_FX_BR.bam \
  -st tumor tumor -ss male male -sg C-52YNHF C-52YNHF -sn C-52YNHF-N001-d C-52YNHF-N002-d \
  --vcf /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/fingerprinting_snps.vcf \
  -db /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/db/ \
  --fafile /Users/charlesmurphy/Desktop/data/ref/hg19/Homo_sapiens_assembly19.fasta \
  --outdir ~/Desktop --json --plot

biometrics major \
  -sb /Users/charlesmurphy/Desktop/mskcc-analyses/200608_compare_qc_tools/manually_count_bases/C-52YNHF-N001-d_cl_aln_srt_MD_IR_FX_BR.bam \
  -st tumor -ss male -sg C-52YNHF -sn C-52YNHF-N001-d \
  --vcf /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/fingerprinting_snps.vcf \
  -db /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/db/ \
  --fafile /Users/charlesmurphy/Desktop/data/ref/hg19/Homo_sapiens_assembly19.fasta \
  --outdir ~/Desktop --json

biometrics minor \
  -sb /Users/charlesmurphy/Desktop/mskcc-analyses/200608_compare_qc_tools/manually_count_bases/C-52YNHF-N001-d_cl_aln_srt_MD_IR_FX_BR.bam \
  -st tumor -ss male -sg C-52YNHF -sn C-52YNHF-N001-d \
  --vcf /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/fingerprinting_snps.vcf \
  -db /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/db/ \
  --fafile /Users/charlesmurphy/Desktop/data/ref/hg19/Homo_sapiens_assembly19.fasta \
  --outdir ~/Desktop --json

biometrics sexmismatch \
  -sb /Users/charlesmurphy/Desktop/mskcc-analyses/200608_compare_qc_tools/manually_count_bases/C-52YNHF-N001-d_cl_aln_srt_MD_IR_FX_BR.bam \
  -st tumor -ss male -sg C-52YNHF -sn C-52YNHF-N001-d \
  --vcf /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/fingerprinting_snps.vcf \
  --bed /Users/charlesmurphy/Desktop/data/innovation/resources/MSK-ACCESS-v1.0/MSK-ACCESS-v1_0-probe-A.sorted.bed \
  -db /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/db/ \
  --fafile /Users/charlesmurphy/Desktop/data/ref/hg19/Homo_sapiens_assembly19.fasta \
  --outdir ~/Desktop --json


for i in `seq 1 1000`
do
  python ../../biometrics/cli.py extract \
    -sb /Users/charlesmurphy/Desktop/mskcc-analyses/200608_compare_qc_tools/manually_count_bases/C-52YNHF-N001-d_cl_aln_srt_MD_IR_FX_BR.bam \
    -st tumor -ss male -sg C-52YNHF -sn C-52YNHF-N002-d-$i \
    --vcf /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/fingerprinting_snps.vcf \
    -db /Users/charlesmurphy/Desktop/mskcc-analyses/201013_fingerprinting/db3/ \
    --fafile /Users/charlesmurphy/Desktop/data/ref/hg19/Homo_sapiens_assembly19.fasta
done
