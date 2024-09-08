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

# 3 局部比对 Local alignment

`makeblastdb -in tomato.p_ctg.fasta -dbtype nucl`  

>-in参数指定参考基因组。  
>-dbtype参数指定数据库类型，nucl为核酸，prot为蛋白质。  

`blastn -query query.fa -db tomato.p_ctg.fasta -out query.fa.blastn \`  
`-evalue 1e-5 -num_threads 10 -outfmt 6\`   

>-query：指定待比对的序列文件（FASTA格式）。  
>-db：指定要比对的数据库。  
>-out：输出结果文件名称。  
>-evalue：设置E值阈值，过滤掉E值大于1e-5的比对结果。
>-num_threads：指定程序使用的CPU线程数量。  
>-outfmt 6：设置输出格式为表格格式，常用于解析和分析。  

>Query：查询序列。  
>Subject：数据库中匹配到的序列。  
>Identity：两个序列在比对区域中的相似性百分比。  
>E-value：期望值，表示在数据库中找到的相似性是随机出现的概率，值越小说明结果越显著。  
>Score：比对得分，反映了比对的质量。  

`# Query id   Subject id    % identity  alignment length mismatches gap opens q.start q.end s.start s.end evalue bit score`  
`seq1         gi|12345678|    99.12       100             1         0        1        100    500      600    1e-50    200`  
`seq1         gi|87654321|    95.00       100             5         0        1        100    200      300    1e-20    180`  
`seq2         gi|24681357|    90.50       80              8         0        10       90     700      780    1e-10    150`  

>Query id：查询序列的ID（在FASTA文件中指定）。  
>Subject id：数据库中匹配的序列ID（通常包含基因信息的标识符）。  
>% identity：比对区域中查询序列与目标序列之间的相似性百分比。  
>alignment length：比对的长度，即匹配区域中包含的碱基或氨基酸数目。  
>mismatches：比对过程中不匹配的碱基/氨基酸数目。  
>gap opens：比对中发生插入或缺失的次数（即gap的个数）。  
>q.start 和 q.end：查询序列中比对区域的起始和终止位置。  
>s.start 和 s.end：目标序列（数据库序列）中比对区域的起始和终止位置。  
>evalue：E值，表示这种比对结果随机出现的概率，值越小越显著。  
>bit score：比对的得分，得分越高表示比对质量越好。  

# 4 全局比对 Global alignment

`mafft --thread 52 --auto sequences.fasta > aligned_output.fasta`

>--thread：指定程序使用的CPU线程数量。  
>--auto：自动选择最合适的比对算法，适用于中小规模的序列集。  
>sequences.fasta：输入的 FASTA 格式的序列文件。  
>aligned_output.fasta：输出文件，包含比对后的序列。  

# 5 高通量数据比对 HTS mapping

`bwa index reference.fasta`

`bwa mem -t 52 reference.fasta reads.1.fastq.gz reads.1.fastq.gz \`  
`| samtools sort -O bam -@ 52 > aligned_reads.sort.bam`

>mem：使用BWA-MEM算法。  
>-t：指定程序使用的CPU线程数量。  
>reference.fasta：参考基因组文件。  
>reads.1.fastq.gz reads.1.fastq.gz：双端测序数据的FASTQ文件，使用gzip进行压缩。  
>管道后的samtools sort：便于后续分析（如变异检测），需对BAM文件按染色体位置进行排序。  
>-O bam：输出格式为BAM。  
>-@：指定程序使用的CPU线程数量。  

`samtools index aligned_reads.sort.bam`

`samtools flagstat aligned_reads.sort.bam`


