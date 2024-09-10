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

# 5 高通量测序数据比对 HTS mapping

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

# 6 基因组比对 Pairwise genome alignment

`minimap2 -t 52 -x asm5 SLL_Heinz.fasta tomato.p_ctg.fasta > tomato.p_ctg_Heinz.paf`

>-t 52：使用52个CPU线程进行并行计算。  
>-x asm5：指定全基因组比对模式，适用于序列差异在0.1%左右的两个基因组。  
>SLL_Heinz.fasta：指定参考基因组序列。  
>tomato.p_ctg_Heinz.paf：输出比对结果，格式为PAF。  



`source activate pafCoordsDotPlotly`

>激活pafCoordsDotPlotly环境

`pafCoordsDotPlotly.R -i tomato.p_ctg_Heinz.paf -o tomato.p_ctg_Heinz_dotplot \`  
`-s -m 10000 -q 10000 -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12`  

>-i：指定上一步PAF格式的比对结果。  
>-o：指定输出结果的文件前缀。  
>-s：根据序列相似性指定点图的颜色。  
>-m：仅显示比对长度大于10000bp的结果。  
>-q：仅显示比对长度大于10000bp的query序列。  
>-r：指定横轴要显示的参考基因组上的染色体。  

# 7 遗传变异的鉴定 Variant calling

`bwa index SLL_Heinz.fasta`

>对SLL_Heinz.fasta番茄参考基因组构建索引。

`bwa mem -t 52 SLL_Heinz.fasta tomato.1_clean.fq.gz tomato.2_clean.fq.gz | \`  
`samtools sort -O bam -@ 52 > tomato.sort.bam`

>-t 52：使用52个CPU线程进行并行计算。  
>tomato.1_clean.fq.gz tomato.2_clean.fq.gz：需要鉴定变异的番茄样品，NGS双端测序数据。  
>-O bam：输出格式为BAM。  
>-@ 52：使用52个CPU线程进行并行计算。


`gatk MarkDuplicates -I tomato.sort.bam -O tomato.sort.markdup.bam \`  
`-M tomato.sort.markdup.metrics.txt`

>-I：需要标记重复的BAM格式比对结果。  
>-O：标记重复后的BAM。  
>-M：检测到的重复条目。  

`samtools index tomato.sort.markdup.bam`  

`samtools faidx SLL_Heinz.fasta`
>构建参考基因组的FASTA索引。

`bcftools mpileup -a FORMAT/SP --threads 52 -Ou -f SLL_Heinz.fasta \`  
`tomato.sort.markdup.bam | bcftools call -m -v > tomato.bcftools.vcf`

>-a FORMAT/SP：在VCF结果文件中FORMAT一列输出SP\
>（Phred-scaled strand bias P-value）信息，用于后续变异的过滤。  
>--threads 52：使用52个CPU线程进行并行计算。  
>-Ou：指定输出格式为非压缩的BCF。  
>-f：指定参考基因组，需要FASTA格式索引。  
>-m：使用multiallelic and rare-variant calling的算法进行变异鉴定。  
>-v：仅输出变异位点。  

```
DP=`samtools depth -a tomato.sort.markdup.bam | awk '{i+=$3}END{print i/NR}'`
```

>计算每个基因组位置上的比对覆盖度，输出基因组有效的平均测序深度，用于变异过滤。  
>-a：计算每个位置的覆盖度，包括覆盖度为0的位点。

`bcftools view --threads 52 -m2 -M2 -e "QUAL < 10 || MQ < 10 || DP < 3 || \`  
`DP>3*$DP || MQBZ < -(3.5+4*DP/QUAL) || RPBZ > (3+3*DP/QUAL) || \`  
`RPBZ < -(3+3*DP/QUAL) || FORMAT/SP > (40+DP/2) || SCBZ > (2.5+DP/30)" \`  
`tomato.bcftools.vcf | bcftools norm --threads 52 -f SLL_Heinz.fasta \`  
`-c x -D > tomato.bcftools.filter.vcf`  

>bcftools view：用于查看、过滤或转换 VCF 或 BCF 格式的文件。  
--threads 52：使用52个CPU线程来并行处理数据。  
-m2 -M2：仅保留二等位变异，不考虑多等位变异。  
-e：表示用表达式对变异进行过滤。如果表达式为真，该变异将被过滤掉。  

>详细解释过滤条件：  
QUAL < 10：过滤质量分数小于 10 的变异，质量分数越高，变异越可信。  
MQ < 10：过滤测序的比对质量（Mapping Quality）小于 10 的变异。  
DP < 3：过滤测序深度（Depth）小于 3 的变异。  
DP > 3*$DP：过滤测序深度超过 3 倍平均有效深度（$DP，上一步计算得出）的变异。  
MQBZ < -(3.5+4*DP/QUAL)：过滤 Mapping Quality Z-score (MQBZ)   
的值低于动态阈值的变异，MQBZ 用于检测比对质量分数的偏差。  
RPBZ > (3+3*DP/QUAL)：过滤 Read Position Bias Z-score (RPBZ)   
的值大于动态阈值的变异，RPBZ 检测测序读长在比对时发生偏位的可能性。  
RPBZ < -(3+3*DP/QUAL)：过滤 RPBZ 值小于该阈值的变异，同样检测测序偏差。  
FORMAT/SP > (40+DP/2)：过滤 Strand Bias（测序方向偏差）大于阈值的变异，  
这种偏差可能表明变异是由测序错误引起的。  
SCBZ > (2.5+DP/30)：过滤 Soft-Clip Bias Z-score (SCBZ) 的值大于阈值的变异，  
SCBZ 用于检测读长是否有被“软裁剪”现象。  

>bcftools norm：这个命令用于对 VCF 文件进行标准化处理（normalization），  
>确保变异格式统一，以便后续分析。  
--threads 52：同样使用52个CPU线程进行并行化处理。  
-f：提供参考基因组文件，用于校正变异位点。  
-c x：去除REF的基因型与参考基因组不一致的位点。  
-D：处理变异中的重复位点，确保同一基因组位置不会被多次记录。  

# 8 系统发生树的构建 Phylogenetic tree

`mafft --auto protein_sequences.fasta > protein_sequences_aligned.fasta`

>--auto：自动选择最适合的比对策略。  
protein_sequences.fasta：输入的蛋白序列文件（FASTA格式）。  
protein_sequences_aligned.fasta：输出的比对结果文件。  

`raxmlHPC -f a -s protein_sequences_aligned.fasta -n tree -m PROTGAMMAAUTO \`  
`-x 12345 -p 12345 -# 100 -T 20`

>-f a：指定分析类型。a 表示进行全面的分析，包括最大似然树的构建和 Bootstrap 置信度评估。  
-s protein_sequences_aligned.fasta：输入的比对文件。  
-n tree：输出文件的前缀。  
-m PROTGAMMAAUTO：选择蛋白质进化模型，PROTGAMMAAUTO 代表自动选择适合的进化模型。  
-x 12345：指定随机种子，用于在Bootstrap分析中生成随机样本的起始点。  
-p 12345：指定随机种子，用于启动树重建过程中的随机数生成。  
-# 100：执行100次Bootstrap重采样，以评估树的可靠性。  
-T 20：指定使用的线程数。将使用20个CPU线程来加速分析过程。


>RAxML_bestTree.tree：最佳的最大似然树。  
RAxML_bootstrap.tree：Bootstrap支持度树。

>可以使用树可视化工具（如FigTree或iTOL）来查看和解释构建的系统发生树。




