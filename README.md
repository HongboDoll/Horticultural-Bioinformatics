# Horticultural-Bioinformatics
Codes used in the book of horticultural bioinformatics

# 1 基因组拼接 Genome assembly

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

# 2 基因组注释 Genome annotation

`EDTA.pl --genome tomato.p_ctg.fasta --step all --overwrite 1 --anno 1 --threads 52`

>--genome参数指定需要注释的参考基因组。--step all参数表示运行全部步骤。  
>--overwrite 1参数表示覆盖上次运行的结果。--anno 1参数会生成一系列后续需要用到的结果文件。  
>--threads参数指定程序使用的CPU线程数量。  

`bedtools maskfasta -soft -fi tomato.p_ctg.fasta \`  
`-fo tomato.p_ctg.fasta.softmasked -bed tomato.p_ctg.fasta.mod.EDTA.TEanno.gff3`  

>使用BEDTools生成soft-masked（软屏蔽重复）的参考基因组。用于后续编码基因预测。

>*EDTA.TEanno.sum文件为注释出重复序列的统计。结果如下：

`hisat2-build -p 52 tomato.p_ctg.fasta tomato.p_ctg.fasta`

>对参考基因组构建索引用于下一步比对。  
>-p参数指定程序使用的CPU线程数量。  

`hisat2 -x tomato.p_ctg.fasta -1 tomato.1_clean.fq.gz -2 tomato.2_clean.fq.gz \`  
`-p 52 | samtools sort - > tomato.sort.bam`

>将转录组双端测序得到的序列比对到参考基因组上。-x参数指定上一步建好的索引。  
>-1和-2参数分别指定双端测序文件。-p参数指定程序使用的CPU线程数量。  
>使用samtools对比对结果进行排序。  

`braker.pl --genome=tomato.p_ctg.fasta.softmasked --workingdir=tomato_braker \`  
`--gff3 --threads=8 --species=tomato --prot_seq=homolog_protein.fa \`  
`--bam=tomato.sort.bam`  

>--genome参数指定软屏蔽重复的参考基因组。--workingdir参数指定程序输出目录。  
>指定--gff3参数会生成GFF3格式的结果。--threads参数指定程序使用的CPU线程数量。  
>--species参数指定生成的AUGUSTUS训练集的名称。--prot_seq参数指定同源蛋白序列。  
>--bam参数指定转录组数据的比对结果。

最后的结果文件为`tomato_braker/braker.gff3`


