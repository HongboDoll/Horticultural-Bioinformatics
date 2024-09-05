# Horticultural-Bioinformatics
Codes used in the book of horticultural bioinformatics

# 1 Genome assembly 基因组拼接

`samtools fastq tomato.hifi_reads.bam | gzip > tomato.ccs.fastq.gz`

`hifiasm --primary -t 52 tomato.ccs.fastq.gz -o tomato`

>当使用--primary参数时，hifiasm只会输出一套代表性的主基因组（primary assembly）  
>而不会输出两套单倍型序列。  
>在基因组杂合性较高的物种中，该参数可以避免生成多套单倍型组装，简化分析。

>-t参数指定程序使用的CPU线程数量。

`awk '/^S/{print ">"$2;print $3}' tomato.p_ctg.gfa > tomato.p_ctg.fasta`

`source activate merqury`

>激活merqury的conda环境

`$MERQURY/best_k.sh 800000000 0.001`

>在番茄估计基因组大小为800Mbp时，计算下一步所需的最佳k值

`meryl k=19 threads=52 count tomato.ccs.fastq.gz output tomato.meryl`

>当k值为19时，使用52个CPU线程构建PacBio原始测序数据的meryl库

`merqury.sh tomato.meryl tomato.p_ctg.fasta tomato.merqury > tomato.merqury.log`

>对拼接质量进行评估
