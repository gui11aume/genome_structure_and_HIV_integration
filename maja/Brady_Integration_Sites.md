PATHTOFastxtoolkit0013=/YOURPATH/fastx_toolkit_0.0.13  
PATHTObbmap=/YOURPATH/bbmap  
PATHTOpblat=/YOURPATH/icebert-pblat-e05e284/pblat  


##### Trimming:
```{bash}

$PATHTObbmap/bbduk2.sh in=SRR8860863.fastq out=SRR8860863_clean.fastq lliteral=TGCATCGAAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA minlength=0 \
threads=12 rliteral=TAGTCCCTTAAGCGGAGCCCTATAGTGCTGATGGCGCGAGGGAGGC editdistance=2

$PATHTObbmap/bbduk2.sh in=SRR8860864.fastq out=SRR8860864_clean.fastq lliteral=TGCATACGAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA minlength=0 \
threads=12 rliteral=TAGTCCCTTAAGCGGAGCCCTATAGTGCTGATGGCGCGAGGGAGGC editdistance=2  

$PATHTObbmap/bbduk2.sh in=SRR8860865.fastq out=SRR8860865_clean.fastq lliteral=TGCAGTCAAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA minlength=0 \
threads=12 rliteral=TAGTCCCTTAAGCGGAGCCCTATAGTGCTGATGGCGCGAGGGAGGC editdistance=2  

$PATHTObbmap/bbduk2.sh in=SRR8860866.fastq out=SRR8860866_clean.fastq lliteral=ATCGCGATAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA minlength=0 \
threads=12 rliteral=TAGTCCCTTAAGCGGAGCCCTATAGTGCTGATGGCGCGAGGGAGGC editdistance=2  
```

#### Fastq -> Fasta 

```{bash}
for i in 3 4 5 6
do
$PATHTOFastxtoolkit0013/bin/fastq_to_fasta -i SRR886086${i}_clean.fastq -o SRR886086${i}_clean.fasta -n -Q33
done
```

##### BLAT to hg19:

```{bash}


for i in 3 4 5 6
do 
$PATHTOpblat -t=dna -q=dna -threads=10 -maxIntron=0 -minIdentity=98 hg19.fa SRR886086${i}_clean.fasta SRR886086${i}onhg19_98.psl
done
```
##### R code:
```{r}
library(data.table)


doItForgenome <- function(genome){
	nms <- c("Resting", "RestingdNTP", "Activated", "ActivateddNTP")

	for (i in 1:4) {
	assign(nms[i],{
		x <- fread(paste(c("SRR886086",i+2,"onhg",genome,"_98.psl"),collapse=""),skip=5)
		x<- x[,Pid:=V1/V11][Pid>=0.98][,isUniquelyMapped:={.N==1},V10]
		x[isUniquelyMapped==TRUE][,.(V10,V14,V16)][,.N,.(V14,V16)]
		})
	}


	afterBlat98 <-sapply(c("Resting", "RestingdNTP", "Activated", "ActivateddNTP"), function(y){
		x <- get(y)
		c(TotalUniquelyMappableIS=x[,sum(N)],TotalUniqueIS=x[,.N])
	})
	t(afterBlat98)

}

doItForgenome(19)

```
#### Results:
```
> doItForgenome(19)
              TotalUniquelyMappableIS TotalUniqueIS
Resting                           963           701
RestingdNTP                       478           386
Activated                         594           363
ActivateddNTP                     881           501




```
