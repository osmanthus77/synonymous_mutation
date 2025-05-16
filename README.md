# 同义/非同义突变计算

- 定义：
  > 同义突变率 ds/Ks：发生同义突变的位点数目（Ds）与可能发生同义突变的所有位点总数的比值
  >
  > 非同义突变率 dN/Ka：发生非同义突变的位点数目（Dn）与可能发生非同义突变的所有位点总数的比值
  >
  > Ka/Ks:非同义替换没有改变蛋白质的组成，因此不受自然选择的影响(忽略密码子偏好性)，那么 Ks 就能反映进化过程的背景碱基替换率。Ka/Ks 的比值就能说明这个基因是受到了何种选择。

- 选择压力
  > 正选择：Ka/Ks>>1，基因受到正选择，在快速进化中得以保留，造成进化上的优势
  >
  > 中性选择：Ka/Ks=1，基因中性进化
  >
  > 负选择/纯化选择：Ka/Ks<<1，基因受到负选择，选择消除有害突变

参考：https://blog.csdn.net/qq_36608036/article/details/131557801

## 0 软件安装

### 0.1 antismash v7.1.0

- 安装`conda`

```shell
# 下载安装conda
mkdir ~/miniconda3
# 下载链接的安装包为intel Mac版本，可根据自己电脑不同下载其他版本
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
source ~/miniconda3/bin/activate

# 初始化并配置环境
conda init --all

# 配置软件源channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

- 使用`conda`手动安装 antismash
参考官方网页 [下载教程](https://docs.antismash.secondarymetabolites.org/install/)

```shell
# 创建环境、激活
conda create -n antismash
conda activate antismash
# 下载依赖
conda install hmmer2 hmmer diamond fasttree prodigal blast
conda install meme==4.11.2
# 下载安装包、解压、安装
cd /Users/yakamozchan/miniconda3/envs/antismash
wget https://dl.secondarymetabolites.org/releases/7.1.0/antismash-7.1.0.tar.gz
tar -zxf antismash-7.1.0.tar.gz
pip install ./antismash-7.1.0
# 下载数据库
download-antismash-databases
conda deactivate
# 测试安装是否成功
antismash --prepare-data

# 输出下面内容则成功安装
# All prerequisites satisfied

# 运行、用法
source ~/miniconda3/bin/activate
conda activate antismash
antismash my_input.gbk
```

### 0.2 序列比对工具

使用`Homebrew`安装。[homebrew 官网及安装](https://brew.sh)

#### mafft v7.526

```shell
brew install mafft
```

#### muscle v5.1.0

```shell
brew install muscle
```

ps:mafft 比 muscle 运行速度慢，但准确性更高

### 0.3 比对结果修剪——trimal v1.4.1

```shell
brew install triaml
```

### 0.4 构树工具

#### fasttree v2.1.11

```shell
brew install fasttree
```

#### iqtree v1.6.12

```shell
brew install iqtree
```

ps: iqtree 比 fasttree 运行慢，但更加准确

### 0.5 newick-utils

对系统发育树进行处理和可视化的工具

```shell
brew install newick-utils
```

### 0.6 树的可视化

table2itol 脚本，[参考链接](https://www.jianshu.com/p/684b83a5e844)

```shell
git clone https://github.com/mgoeker/table2itol.git
chmod +x table2itol.R
./table2itol.R
```

### 0.7 KaKs_Calculator v3.0

[下载地址](https://github.com/Chenglin20170390/KaKs_Calculator-3.0.git)

```shell
cd ~/biosoft
git clone https://github.com/Chenglin20170390/KaKs_Calculator-3.0.git
cd KaKs_Calculator-3.0/bin
make
echo 'export PATH=$HOME/biosoft/KaKs_Calculator-3.0/bin:$PATH'  >> $HOME/.bashrc
source ~/.bashrc
```

### 0.8 ParaAT v2.0

整合了计算 KaKs 所需的一整套分析，包括蛋白序列比对、根据比对结果回译成 codon 对应的核酸比对结果、计算 KaKs 值（使用 KaKs_Calculator）  
[下载地址 1](https://ngdc.cncb.ac.cn/tools/paraat)  
[下载地址 2](https://github.com/jdebarry/paraat.git)

```shell
# 解压
unzip ParaAT2.0.tar.gz
# 设置环境变量，后续可直接作为命令使用
echo 'export PATH=$HOME/biosoft/ParaAT2.0:$PATH'  >> $HOME/.bashrc
source ~/.bashrc

# 依赖工具下载（蛋白比对）
brew install mafft muscle

# 测试是否安装成功
ParaAT.pl -h test.homologs -n test.cds -a test.pep -p proc -m mafft -f axt -g -k -o result_dir
# 成功安装后有输出目录result_dir
```

#### ParaAT2.0 用法

```shell
ParaAT.pl -h test.homologs -n test.cds -a test.pep -p proc -m mafft -f axt -g -k -o result_dir
```

- `-h`：同源基因名称文件
- `-n`：指定核酸序列文件
- `-a`：指定蛋白序列文件
- `-p`：线程数
- `-m`：指定比对工具（clusterw2、t_coffee、mafft、muscle 其中一个）
- `-g`：去除比对有 gap 的密码子
- `-k`：用 KaKs_Calculator 计算 KaKs 值
- `-o`：输出结果的目录
- `-f`：输出比对文件的格式

### R 包

```R
# pheatmap v1.0.12
# readr v2.1.5
# dplyr v1.1.4

install.packages(c("pheatmap", "readr", "dplyr"))
```

### 0.9 nwr

使用rust的包管理器`cargo`从github安装nwr：

```shell
# 先安装rust及cargo
curl https://sh.rustup.rs -sSf | sh

# 安装nwr
cargo install --git https://github.com/wang-q/nwr
```

### 0.10 tectonic

`tectonic`用于将`.tex`转换为`.psf`，使用`homebrew`进行安装

```shell
brew install tectonic
```

## 1 下载数据

根据 https://github.com/wang-q/genomes/blob/main/groups/Bacillus.md#species-with-assemblies 下载得到类芽孢杆菌 Paenibacillaceae 基因组序列(genbank 文件)，分为有 species 分类和无 species 分类，分别为 complete_taxon、complete_untaxon

## 2 运行 antismash 预测

### 2.1 使用 antismash 注释预测

```shell
mkdir ~/chenxy/Pae_rerun
cd ~/chenxy/Pae_rerun
for level in complete_taxon complete_untaxon; do
    mkdir -p antismas_out/bgc_7.1/Paenibacillaceae/${level}
    find ${level}/gbk/*.gbk |
        xargs baselevel -s .gbk |
        parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 "
            echo {}; antismash \
                --taxon bacteria -c 12 \
                --cb-general --cc-mibig \
                --cb-knownclusters \
                --pfam2go ${level}/gbk/{}.gbk \
                --output-dir antismash_out/bgc_7.1/${level}/{}
            "
```

得到结果目录中`region.html`文件：
![antismash_html_picture](/pic/antismash_html.png "antismash_html")

## 3 提取整理 antismash 预测结果

从`region.html`中`overview`表格中提取预测结果，并从中提取产物信息，再筛选目标产物 polymyxin

### 3.1 统计 strain

先统计使用 antismash 预测的菌株名称，得到 tsv 文件

```shell
cd ~/chenxy/Pae_rerun/result
mkdir strains_raw
for i in Paenibacillaceae; do
    for n in complete_taxon complete_untaxon; do
        find ../bgc_7.1/${i}/${n}  -mindepth 1 -maxdepth 1 -type d  |
            xargs -I {} baselevel {} > strains_raw/strains_${i}_${n}_7.1.lst
        wc -l strains_raw/strains_${i}_${n}_7.1.lst
    done
done
# 输出
# 178 strans_raw/strains_Paenibacillaceae_complete_taxon_7_1.lst
# 154 strans_raw/strains_Paenibacillaceae_complete_untaxon_7_1.lst
```

### 3.2 提取 antismash 结果 html 中产物信息

#### 3.2.1 基于 overview-table 得到所有产物信息

从`region.html`文件中提取`overview`表格的 html 信息，得到每个菌株的原始文件`raw.tsv`，并修改其格式（调整到与 html 界面中 overview table 一致）。最终得到所有菌株的全部产物信息

```shell
# 得到overview的raw文件
cd ~/chenxy/Pae_rerun/result
for level in complete_taxon complete_untaxon ; do
    for version in 7.1; do
        raw_dir=product/product_raw/${level}
        mkdir -p ${raw_dir}

        for strains in $(cat strains_raw/strains_Paenibacillaceae_${level}_${version}.lst); do
            cat ../antismash_out/bgc_${version}/Paenibacillaceae/${level}/${strains}/index.html |
                pup 'table.region-table tbody tr td text{}' |
                sed 's/Region/|Region/g' |
                grep '\S' > ${raw_dir}/${strains}_product_raw.tsv
        done
    done
done

# 修改raw文件格式
mkdir -p  product/product_whole
for level in complete_taxon complete_untaxon ; do
	for version in 7.1; do
		for strains in $(cat strains_raw/strains_Paenibacillaceae_${level}_${version}.lst); do
            raw_dir=product/product_raw/${level}
            perl ../script/html.pl ${raw_dir}/${strains}_product_raw.tsv |
                sed "s/^/${strains}_/g; s/Region /cluster/g" >> product/product_whole/product_whole_bgc_${level}_${version}.tsv
		done
	done
done


## ${raw_dir}/${strains}\_product_raw.tsv 文件结构：

    |Region 1 # Region
    lassopeptide # Type
    694,701 # From
    718,615 # To
    paeninodin # Most similar known cluster
    RiPP # Most similar known cluster
    100% # Similarity

## 结果文件格式
# Paenib_polym_M1_GCF_000237325_1_cluster1	NRPS	67627	131356	fusaricidinB	Polyketide+NRP:Lipopeptide	100%
# Paenib_polym_M1_GCF_000237325_1_cluster2	RRE-containing	263117	283380
# Paenib_polym_M1_GCF_000237325_1_cluster3	transAT-PKSNRPS	1068969	1146361
```

#### 3.2.2 提取 overview 表格中产物的 mibig 信息

从`region.html`文件中提取`overview`表格中提取预测产物的 mibig 参考信息，即 MIBiG 数据库中 BGC 编号

```shell
cd ~/chenxy/Pae_rerun/result
mkdir -p product/product_mibig
for level in complete_taxon complete_untaxon; do
    for version in 7.1; do
        for strains in $(cat strains_raw/strains_Paenibacillaceae_${level}_${version}.lst); do
            cat ../antismash_out/bgc_7.1/Paenibacillaceae/${level}/${strains}/index.html |
                pup 'table.region-table tbody tr.linked-row a attr{href}' |
                grep -v "https://docs.antismash.secondarymetabolites.org/" |
                grep -B 1 "https://mibig.secondarymetabolites.org/" |
                grep -v "\-\-" |
                sed -E 's/https:\/\/mibig.secondarymetabolites.org\/go\///g; s/#//g; s/\/[0-9]//g' |
                paste - - |
                sed "s/^/${strains}_/g" |
                sed -E 's/r1c/cluster/g; s/r([0-9]+)c([0-9]+)/cluster\1.\2/g'
        done > product/product_mibig/mibig_whole_bgc_${level}_${version}.tsv

        # 统计预测产物MiBIG参考信息的个数
        wc -l product/product_mibig/mibig_whole_bgc_${level}_${version}.tsv
    done
done

# 1055 product/product_mibig/mibig_whole_bgc_complete_taxon_7.1.tsv
# 719 product/product_mibig/mibig_whole_bgc_complete_untaxon_7.1.tsv
# 结果格式：Paenib_polym_M1_strain_000237325_1_cluster1	BGC0001152
```

- 统计 mibig 数据库中 polymyxin 参考序列

下载 MIBiG 数据库序列，并统计其中有关目标产物 polymyxin 的序列

```shell
cd ~/chenxy/Pae_rerun/result
mkdir -p mibig/mibig
wget https://dl.secondarymetabolites.org/mibig/mibig_json_4.0.tar.gz -O mibig/MIBIG4_4.0.tar.gz
tar -xvf mibig/MIBIG4_4.0.tar.gz -C mibig/mibig  --strip-components 1
ls mibig/mibig | cut -d "." -f 1 > mibig/mibig.tsv
wc -l mibig/mibig.tsv   #3013

# 统计其中有关polymyxin的序列
cat mibig/mibig.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 '
        cat mibig/mibig/{}.json | 
        grep -i -E "polymyxin|colistin|macolacin" | 
        sed "s/^/{}_/g"
    '

## BGC0000408/BGC0001153/BGC0001192/BGC0002653
## BGC0000408_                "compound": "polymyxin",
## BGC0001153_                    "d-dab polymyxin b"
## BGC0001153_                "compound": "polymyxin B",
## BGC0001192_                "compound": "colistin A",
## BGC0001192_                "compound": "colistin B",
## BGC0002653_                "compound": "macolacin"
```

#### 3.2.3 基于 overview 和 mibig 信息筛选产物 polymyxin

合并 overview 文件和 mibig 文件，筛选预测相似度高于 50%部分，并筛选预测产物为 polymyxin 的。再将筛选到的菌株的 antismash 结果拷贝过来

```shell
# 合并所有对应的 BGCs的预测产物overview信息和 MiBIG参考信息
cd ~/chenxy/Pae_rerun/result
for level in complete_taxon complete_untaxon; do
    for version in 7.1; do
        cat product/product_whole/product_whole_bgc_${level}_${version}.tsv |
            tsv-join -d 1 -f product/product_mibig/mibig_whole_bgc_${level}_${version}.tsv -k 1 --append-fields 2 | # 合并product和MiBIG
            sed 's/%//g' |
            tsv-filter --ge 7:50 | # 筛选：相似度≥50%,可选
            sed 's/_cluster/\tcluster/g' |
            grep -E "polymyxin|colistin|macolacin" > product/polymyxin_${level}_${version}.tsv # 筛选：polymyxin
        wc -l product/polymyxin_${level}_${version}.tsv

    done
done

# 结果文件格式
# Paenib_polym_M1_GCF_000237325_1	cluster16	NRPS	5046972	5128182	polymyxin	NRP	100	BGC0000408
# Paenib_polym_YT9_GCA_039719015_1	cluster15	NRPS	4961917	5043091	polymyxin	NRP	100	BGC0000408
# Paenib_polym_ATCC_842_GCF_022811565_1	cluster8	NRPS	2793810	2873397	polymyxin	NRP	100	BGC0000408

# 拷贝相关antismash结果到新文件夹
for level in complete_taxon complete_untaxon; do
    for version in 7.1; do
        cat product/polymyxin_${level}_${version}.tsv |
            cut -f 1 |
            parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 8 "
                echo {};
                mkdir -p product//polymyxin_${level}_${version}/{};
                cp -r ../antismash_out/bgc_${version}/Paenibacillaceae/${level}/{} product/polymyxin_${level}_${version}/
                "
    done
done
```

从筛选到 polymyxin 的结果文件中选择部分菌株作为**example**进行后续分析：

```tsv
example_polymyxin_complete_taxon_7.1.tsv文件内菌株

Paenib_barc_KACC11450_GCF_013347305_1	cluster5
Paenib_dendri_2022CK_00834_GCF_029952845_1	cluster10
Paenib_dendri_J27TS7_GCF_021654795_1	cluster16
Paenib_lentus_DSM_25539_GCF_003931855_1	cluster2
Paenib_peo_ZBSF16_GCF_022531965_1	cluster14
Paenib_peo_ZF390_GCF_014692735_1	cluster14
Paenib_polym_CPL258_GCA_029543045_1	cluster2
Paenib_polym_R_5_31_GCF_023586705_2	cluster11
Paenib_polym_R_6_14_GCF_023586785_2	cluster16
Paenib_polym_Sb3_1_GCF_000819665_1	cluster10
Paenib_thi_Mbale2_GCF_028554415_1	cluster15
Paenib_thi_SY20_GCF_028226775_1	cluster17
Paenib_xylane_PAMC_22703_GCF_001908275_1	cluster16
```

## 4 从 antismash 结果提取 domain 序列

提取思路：antismash 结果目录中，`region.js`存储了相关信息。详见[antismash 解读](https://github.com/osmanthus77/something/blob/main/antismash.md)

在`region.js`文件中，从`details_data`中解析`nrpspks`数据，遍历每个`orfs`和`domains`，从中提取出 domain 的 polymyxin 对应的菌株 bgc 的全部 domain 信息

`details_data`部分内容如图：
![regions_js_picture](/pic/regions_js.png "regions_js")

### 4.1 提取 domain DNA 序列

```shell
# 修改第二列格式
cd ~/chenxy/Pae_rerun/result
awk -F'\t' 'BEGIN{OFS="\t"} {gsub("cluster", "r1c", $2); print}' product/example_polymyxin_complete_taxon_7.1.tsv > temp.tsv && mv temp.tsv product/example_polymyxin_complete_taxon_7.1.tsv
cut -f 1,2 product/example_polymyxin_complete_taxon_7.1.tsv > product/example_polymyxin.tsv
head product/example_polymyxin_complete_taxon_7.1.tsv
# 格式
# Paenib_barc_KACC11450_GCF_013347305_1	r1c5
# Paenib_dendri_2022CK_00834_GCF_029952845_1	r1c10
# Paenib_dendri_J27TS7_GCF_021654795_1	r1c16

# 提取全部 domain
mkdir -p domain_dna
for i in $(cat product/example_polymyxin.tsv | sed "s/\t/,/g"); do
    echo ${i};
    sample=$(echo ${i} | cut -d "," -f 1);
    num=$(echo ${i} | cut -d "," -f 2);
    echo ${num};
    js="product/polymyxin_complete_taxon_7.1/${sample}/regions.js";
    type=$(python ../script/antismash_pp_for_complete_uncomplete.py "${js}" "${num}" "${sample}");
    echo ${type} | sed "s/]/]\n/g" >> domain_dna/domain_dna_all_polymyxin_complete_plasmid_panenbacillus.txt;
done
```

提取 domain 后`txt`文件修改格式并筛选特定的 domain，以 C 域为例；  
！注意：筛选特定 domain 时，要用 domain 全称，如下表
| domain 简称 | domain 全称 |
|:----------:|--------------|
| A domain | AMP-binding |
| C domain | Condensation |
| T domain | PCP |

```shell
cd ~/chenxy/Pae_rerun/result
cat domain_dna/domain_dna_all_polymyxin_complete_plasmid_panenbacillus.txt |
sed "s/\[//g" | sed "s/\]//g"|sed "s/'//g" | sed "s/,/\n/g"|
grep -E "Condensation" |       # 筛选特定domain
sed "s/ >//g" | sed "s/+/\t/g" |
sed "s/^/>/g" |
sed "s/\t/\n/g" > domain_dna/domain_Cdna.txt
# 每个菌株第一行 >> 替换为 >
sed "s/^>>/>/g" domain_dna/domain_Cdna.txt > temp.txt && mv temp.txt domain_dna/domain_Cdna.txt
sed "s/^>>/>/g" domain_dna/domain_Adna.txt > temp.txt && mv temp.txt domain_dna/domain_Adna.txt
sed "s/^>>/>/g" domain_dna/domain_Tdna.txt > temp.txt && mv temp.txt domain_dna/domain_Tdna.txt
```

### 4.2 提取 domain 氨基酸序列

方法同`4.1`，将`antismash_pp_for_complete_uncomplete.py`脚本中 dna_sequence 替换为 sequence，即为氨基酸序列。

```shell
# 提取domain 氨基酸序列
mkdir -p domain_aa
for i in $(cat product/example_polymyxin.tsv | sed "s/\t/,/g"); do
    echo ${i};
    sample=$(echo ${i} | cut -d "," -f 1);
    num=$(echo ${i} | cut -d "," -f 2);
    echo ${num};
    js="product/polymyxin_complete_taxon_7.1/${sample}/regions.js";
    type=$(python ../script/antismash_aa_for_complete_uncomplete.py "${js}" "${num}" "${sample}");
    echo ${type} | sed "s/]/]\n/g" >> domain_aa/domain_aa_all_polymyxin_complete_plasmid_panenbacillus.txt;
done

# 修改格式并筛选特定 domain，以C域为例
cat domain_aa/domain_aa_all_polymyxin_complete_plasmid_panenbacillus.txt |
sed "s/\[//g" | sed "s/\]//g"|sed "s/'//g" | sed "s/,/\n/g"|
grep -E "Condensation" |       # 筛选特定domain
sed "s/ >//g" | sed "s/+/\t/g" |
sed "s/^/>/g" |
sed "s/\t/\n/g" > domain_aa/domain_Caa.txt
# 每个菌株第一行 >> 替换为 >
sed "s/^>>/>/g" domain_aa/domain_Caa.txt > temp.txt && mv temp.txt domain_aa/domain_Caa.txt
sed "s/^>>/>/g" domain_aa/domain_Aaa.txt > temp.txt && mv temp.txt domain_aa/domain_Aaa.txt
sed "s/^>>/>/g" domain_aa/domain_Taa.txt > temp.txt && mv temp.txt domain_aa/domain_Taa.txt
```

- 提取 polymyxin 的底物氨基酸

```shell
cd ~/chenxy/Pae_rerun/result
mkdir -p amino_acid_residue/polymyxin
for level in complete_taxon complete_untaxon; do
    for version in 7.1; do
        file_dir="product/product_whole/product_whole_bgc_${level}_${version}.tsv"

        output_file="amino_acid_residue/polymyxin/polymyxin_amino_acid_residue_${level}_${version}.tsv"

        cat "$file_dir"| grep -v '^$' |grep -E "polymyxin|colistin|macolacin"|
        while IFS=$'\t' read -r strain_path_cluster _ _ _ _ _ _ _ _; do
		    strain_path=$(echo ${strain_path_cluster}| sed "s/_cluster/\t/g"| cut -f 1)
            num=$(echo ${strain_path_cluster}| sed "s/_cluster/\t/g"| cut -f 2)
            path="../antismash_out/bgc_${version}/Paenibacillaceae/${level}/${strain_path}/${strain_path}.json"
			html_path="../antismash_out/bgc_${version}/Paenibacillaceae/${level}/${strain_path}/index.html"
			X_aa=$(cat  ${html_path} | pup "div#r1c${num} dl.prediction-text " | grep -v '^$''^$'| grep -E "^:|nrpys:" |sed -E "s/\s+//g" |sed -E "s/nrpys:/(/g" | tr "\n" ")")
            echo -e "${strain_path}\t${num}\t${level}_${version}\t${X_aa}"
        done > "$output_file"
    done
done
```

## 5 比对构树

### 5.1 多序列比对

多序列比对（Multiple sequence alignment, MSA）：将多条有系统进化关系的氨基酸序列 or 核酸序列进行比对，把相同的碱基 or 氨基酸残基排在同一列上——对齐。做序列对齐的主要目的是，确定所有序列的同源位点相互对应。

```shell
cd ~/chenxy/Pae_rerun/result
mkdir tree
mafft --auto domain_dna/domain_Cdna.txt > tree/domain_Cdna.aln.fa
```

`--auto`选项：自动选择比对模式为高速模式 or 高精度模式。高精度适用于少于约 200 个序列 x 每个序列少于约 2000 个氨基酸/核苷酸。

### 5.2 多序列比对结果修剪

目前存在各种各样的多序列比对算法，但是不存在一个算法能够绝对地保证其能进行完美的位点对应。此外，用于比对的序列可能存在一些错误或者删除和缺失。所以，在得到多序列比对结果之后，用于进化树构建之前，对多序列比对结果进行修剪。

- 修剪操作：

  > 人工修剪
  >
  > 删除所有含有 gaps 的位点
  >
  > 删除不保守位点，如 G-blocks 软件
  >
  > 删除含有一定比例的 gaps 的位点，如 MEGA 内置算法
  >
  > 按照位点信息量删除，如 trimAL

- trimal 用法：

```shell
trimal -in <input> -out <output> -[option]
```

常用 option：`-automated1`自动选择最优修剪方式，优化用于最大似然法构树

一般使用`trimAL`进行修剪:

```shell
cd ~/chenxy/Pae_rerun/result
trimal -in tree/domain_Cdna.aln.fa -out tree/domain_Cdna.trim.fa -automated1
```

### 5.3 构树

- 系统发育树分为有根树、无根树
  > 无根树只反应分类单元之间的距离、无祖先问题
  >
  > 有根树反映进化时间顺序，树枝长度可以反应不同基因/蛋白进化方式和进化速率

- 常见算法：
  > 最大似然法(Maximum Likelihood, ML)
  >
  > 邻接法(Neighbor-Joining, NJ)
  >
  > 最大简约法(Maximum Parsimony, MP)
  >
  > 最小进化法(Minimum Evolution, ME)
  >
  > 贝叶斯推断(Bayesian)
  >
  > 非加权分组平均法(UPGMA)

一般用最大似然法即可，贝叶斯法最准确，但运行速度最慢。  
常用工具：iqtree、fasttree，fasttree 适合快速构树，iqtree 准确性更高。  
详见 https://github.com/osmanthus77/something/blob/main/record.md 中 iqtree 和 fasttree 部分。

这里，使用`fasttree`构树

```shell
cd cd ~/chenxy/Pae_rerun/result
mkdir -p tree/fasttree
fasttree -nt tree/domain_Cdna.trim.fa > tree/fasttree/domain_Cdna.nwk
```

注：`-nt`指定核酸序列，默认为蛋白。

### 5.4 对树操作

#### 5.4.1 简化树

合并具有相同拓扑结构的子树，减少树的复杂度

```shell
cd cd ~/chenxy/Pae_rerun/result
nw_condense tree/fasttree/domain_Cdna.nwk > tree/fasttree/domain_Cdna.condense.nwk
```

#### 5.4.2 可视化

先手动生成注释文件，再使用`table2itol`脚本进行可视化

- 生成注释文件

```shell
# 提取发育树节点名称
cd ~/chenxy/Pae_rerun/result
cat domain_dna/domain_Cdna.txt | grep ">" | sed "s/^>//g" > domain_dna/domain_Cdna_strain.txt
head -n 3 domain_dna/domain_Cdna_strain.txt
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_Starter-6-298-1.1-none
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_LCL-1070-1355-1.4-none
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_LCL-2120-2404-1.7-none


# 手动注释
# 对domain_Cdna_strain.txt中 Cdomain 手动筛选对应菌株antismash结果index.html文件中polymyxin基因、C域类型、C域在每个基因中的位置
# 得到文件domain_Cdna_annotate.txt
head -n 3 domain_dna/domain_Cdna_annotate.txt
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_Starter-6-298-1.1-none	pmxE|C_Starter
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_LCL-1070-1355-1.4-none	pmxE|LCL|Location2
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_LCL-2120-2404-1.7-none	pmxE|LCL|Location3


# 再注释C域类型（C_Starter、LCL、DCL）
tsv-select -f 2 domain_dna/domain_Cdna_annotate.txt | cut -d "|" -f 2 > domain_dna/temp.txt
paste domain_dna/domain_Cdna_annotate.txt domain_dna/temp.txt > domain_dna/domain_Cdna_annotation.txt
# 加上haeder
awk 'BEGIN {print "strain\tcircle1\tcircle2"} {print}' domain_dna/domain_Cdna_annotation.txt > domain_dna/domain_Cdna_annotation_header.txt

# domain_Cdna_annotation_header.txt即为后续用到的注释文件
head -n 3 domain_dna/domain_Cdna_annotation_header.txt
# strain	circle1	circle2
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_Starter-6-298-1.1-none	pmxE|C_Starter	C_Starter
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_LCL-1070-1355-1.4-none	pmxE|LCL|Location2	LCL
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_LCL-2120-2404-1.7-none	pmxE|LCL|Location3	LCL
```

- **重新建树**

由于注释过程中，存在部分菌株 bgc 中预测的 polymyxin 的 nrps 基因多于 mibig 库中，在手动注释时删除这部分内容，如下图：
![zhushishanchu](/pic/注释删除.png "zhushishanchu")
图中，`Query`部分为菌株 bgc 与 mibig 数据库中已知 nrps 基因比对，其中有灰色线条连接部分才是菌株 bgc 中有关 polymyxin 的 nprs 基因，但在`4.1`的提取过程中左边三个无关的红色基因同样也提取出来了，因此需要删除这部分 C 域的序列后重新建树。

```shell
# 删除得到正确的Cdomain 序列
cd ~/chenxy/Pae_rerun/result
cut -f 1 domain_dna/domain_Cdna_annotation.txt > domain_dna/tmp.txt
grep -A 1 -f domain_dna/tmp.txt domain_dna/domain_Cdna.txt > domain_dna/domain_Cdna_new.txt
sed -e "s/--//g" -e "/^$/d" domain_dna/domain_Cdna_new.txt > domain_dna/temp.txt && mv domain_dna/temp.txt domain_dna/domain_Cdna_new.txt
rm domain_dna/tmp.txt

faops size domain_dna/domain_Cdna.txt | wc -l
# 132
faops size domain_dna/domain_Cdna_new.txt | wc -l
# 121
# 序列数量减少，确认删除成功

# 重新建树
mafft --auto domain_dna/domain_Cdna_new.txt > tree/domain_Cdna_new.aln.fa
trimal -in tree/domain_Cdna_new.aln.fa -out tree/domain_Cdna_new.trim.fa -automated1
fasttree -nt tree/domain_Cdna_new.trim.fa > tree/fasttree/domain_Cdna_new.nwk
```

- 使用`table2itol`脚本

```shell
cd ~/chenxy/Pae_rerun/result
mkdir -p tree/fasttree/style1
script=../script/table2itol/table2itol.R
Rscript ${script} -a -D tree/fasttree/style1 -i strain -l Location -w 0.5 domain_dna/domain_Cdna_annotation_header.txt
Rscript ${script} -a -D tree/fasttree/style1 -i strain -l Type -w 0.5 domain_dna/domain_Cdna_annotation_header.txt
```

参数:

> `-a`：表示如果找不到输入列，则程序会终止运行（默认情况下不会执行）
>
> `-D`：指定输出目录
>
> `-i`：指定发育树叶片节点的名称，即注释文件中第一列
>
> `-l`：指定注释信息的名称，即注释文件的第二列 or 第三列

最终在`tree/fasttree/style1`目录下得到`circle1的`两个文件和`circle2`的两个文件，文件结构如下：

```shell
tree/fasttree
├── domain_Cdna.condense.nwk
├── domain_Cdna.nwk
└── style1
    ├── iTOL_colorstrip-circle1.txt
    ├── iTOL_colorstrip-circle2.txt
    ├── iTOL_labels-circle1.txt
    └── iTOL_labels-circle2.txt
```

- iTOL 可视化

> 在[iTOL 网站](https://itol.embl.de)注册账号并上传发育树文件`domain_Cdna.condense.nwk`
>
> 打开该树的文件后，将`tree/fasttree/style1`目录中`iTOL_colorstrip-circle1.txt`和`iTOL_colorstrip-circle2.txt`拖进 itol 页面即完成可视化
>
> 最后在`Export`导出

![fasttree](/pic/fasttree.png "fasttree")

## 6 Ka/Ks 计算

计算 Ka/Ks 需要用到`ParaAT2.0`和`KaKs_Calculator`两个工具。其中需要用到以下三个文件：

> `-h`：同源基因名称文件
>
> `-n`：指定核酸序列文件，即`5.4.2`步骤中`重新建树`筛选出的 dna 序列
>
> `-a`：指定蛋白序列文件

ps：在`4.2`中提取到了 Cdomain 氨基酸序列，但是同样地，由于注释原因删除了部分核酸序列，因此之前提取的氨基酸序列也要进行重新筛选

### 6.1 准备需要的核酸和蛋白序列文件

#### 6.1.1 Cdomain 氨基酸序列筛选

```shell
cd ~/chenxy/Pae_rerun/result
cut -f 1 domain_dna/domain_Cdna_annotation.txt > domain_aa/tmp.txt
grep -A 1 -f domain_aa/tmp.txt domain_aa/domain_Caa.txt > domain_aa/domain_Caa_new.txt
sed -e "s/--//g" -e "/^$/d" domain_aa/domain_Caa_new.txt > domain_aa/temp.txt && mv domain_aa/temp.txt domain_aa/domain_Caa_new.txt
rm domain_aa/tmp.txt

faops size domain_aa/domain_Caa.txt | wc -l
# 132
faops size domain_aa/domain_Caa_new.txt | wc -l
# 121
# 序列数量减少，确认删除成功。并且序列数量与对应的dna序列数量一致，均为121
# 成功生成筛选后的蛋白序列文件
```

#### 6.1.2 将 C domain 序列按位置 Location 拆分为 10 文件

Polymyxin 多粘菌素为十个氨基酸组成的多肽，其 NRPS 模块为十个，其菌株 BGCs 中 C domain 分别位于基因 pmxA 和 pmxE 上，一共十个 C domain，也就是有十个 Location。在此，进行 KaKs 需要以下两个维度进行计算：

> 同一个 BGC 内：某一条特定菌株 BGC 中，计算十个 Location 上的 Cdomain 两两之间的 KaKs
>
> 同一个 Location：针对某一个特定 Location，计算所有菌株 BGC 中该 Location 上的 Cdomain 两两之间的 KaKs

因此，需要先根据注释文件中 Cdomai 信息（基因、类型、位置，如`pmxE|LCL|Location2`），以基因和位置为依据，进行拆分，十个位置拆分为十个文件

- 先生成位置文件 lst

```shell
# 先生成位置文件lst
cd  ~/chenxy/Pae_rerun/result
mkdir KaKs
cut -f 2 domain_dna/domain_Cdna_annotation.txt | sort | uniq > KaKs/10_Location.lst
# pmxA|DCL|Location1
# pmxA|LCL|Location2
# pmxA|LCL|Location3
# pmxA|LCL|Location4
# pmxE|C_Starter
# pmxE|DCL|Location4
# pmxE|LCL|Location2
# pmxE|LCL|Location3
# pmxE|LCL|Location4
# pmxE|LCL|Location5
# pmxE|LCL|Location6
# 手动修改10_Location.lst内容，将Location4的Cdomain两种类型DCL和LCL合并，并删除两种类型

# pmxA|DCL|Location1
# pmxA|LCL|Location2
# pmxA|LCL|Location3
# pmxA|LCL|Location4
# pmxE|C_Starter
# pmxE|C|Location4
# pmxE|LCL|Location2
# pmxE|LCL|Location3
# pmxE|LCL|Location5
# pmxE|LCL|Location6
```

- 按照 Location 拆分完整的核酸和蛋白序列文件为 10 个 Location 序列文件

```shell
# 先将注释文件拆分为十个文件
cd  ~/chenxy/Pae_rerun/result
mkdir KaKs/location
for location in $(cat KaKs/10_Location.lst); do
    cat domain_dna/domain_Cdna_annotation.txt |
    grep "${location}" > KaKs/location/${location}.lst
done
# 删除十个文件中的C domain注释信息，只保留第一列名称
for location in $(cat KaKs/10_Location.lst); do
    cut -f 1 KaKs/location/${location}.lst > KaKs/${location}_tmp.txt;
    mv KaKs/${location}_tmp.txt KaKs/location/${location}.lst;
done

# 根据locaotion.lst文件 拆分核酸序列
mkdir KaKs/Cdna_fasta
for location in $(cat KaKs/10_Location.lst); do
    grep -A 1 -f KaKs/location/${location}.lst domain_dna/domain_Cdna_new.txt >> KaKs/Cdna_fasta/${location}_Cdna.fa;
    sed -e "s/--//g" -e "/^$/d" KaKs/Cdna_fasta/${location}_Cdna.fa > KaKs/${location}_temp.txt && mv KaKs/${location}_temp.txt KaKs/Cdna_fasta/${location}_Cdna.fa;
done

# 根据locaotion.lst文件 拆分氨基酸序列
mkdir KaKs/Caa_fasta
for location in $(cat KaKs/10_Location.lst); do
    grep -A 1 -f KaKs/location/${location}.lst domain_aa/domain_Caa_new.txt >> KaKs/Caa_fasta/${location}_Caa.fa;
    sed -e "s/--//g" -e "/^$/d" KaKs/Caa_fasta/${location}_Caa.fa > KaKs/${location}_temp.txt && mv KaKs/${location}_temp.txt KaKs/Caa_fasta/${location}_Caa.fa;
done

# 确认核酸和氨基酸序列数量与lst文件相同
faops size KaKs/domain_Cdna/pmxA\|DCL\|Location1_Cdna.fa | wc -l
# 12
faops size KaKs/domain_Caa/pmxA\|DCL\|Location1_Caa.fa | wc -l
# 12
cat KaKa/location/pmxA\|DCL\|Location1.lst | wc -l
# 12

```

### 6.2 准备同源基因名称文件

```shell
# 生成一对一对的基因名称文件
# mapfile命令须在 bash4.0 以上版本中使用
cd ~/chenxy/Pae_rerun/result/KaKs
mkdir homologs
for location in $(cat 10_Location.lst); do
    mapfile -t domains < "location/${location}.lst"
    > homologs/${location}.homologs
    for i in "${!domains[@]}"; do
        for j in "${!domains[@]}"; do
            if [ "$i" -lt "$j" ] && [ "${domains[$i]}" != "${domains[$j]}" ]; then
            echo "${domains[$i]},${domains[$j]}" | sed "s/,/\t/g">> homologs/${location}.homologs
            fi
        done
    done
done
head -n 5 homologs/pmxA\|DCL\|Location1.homologs
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_DCL-1320-1619-3.5-none	Paenib_dendri_2022CK_00834_GCF_029952845_1-Condensation_DCL-1337-1636-2.5-none
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_DCL-1320-1619-3.5-none	Paenib_dendri_J27TS7_GCF_021654795_1-Condensation_DCL-1322-1621-3.5-none
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_DCL-1320-1619-3.5-none	Paenib_lentus_DSM_25539_GCF_003931855_1-Condensation_DCL-1304-1603-2.5-none
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_DCL-1320-1619-3.5-none	Paenib_peo_ZBSF16_GCF_022531965_1-Condensation_DCL-1335-1634-3.5-none
# Paenib_barc_KACC11450_GCF_013347305_1-Condensation_DCL-1320-1619-3.5-none	Paenib_peo_ZF390_GCF_014692735_1-Condensation_DCL-1335-1634-10.5-none
```

### 6.3 计算 Ka/Ks 值

使用`ParaAT2.0`脚本进行计算，详细用法见`软件安装`中`0.8`部分

```shell
cd ~/chenxy/Pae_rerun/result/KaKs
# 新建线程数量文件
touch proc
echo "6" > proc

# 计算KaKs
for location in $(cat 10_Location.lst); do
    homologs=homologs/${location}.homologs;
    dna=Cdna_fasta/${location}_Cdna.cds;
    aa=Caa_fasta/${location}_Caa.pep;
    ParaAT.pl -h ${homologs} -n ${dna} -a ${aa} -p proc -m mafft -f axt -g -k -o ${location}_result_dir
done
```

参数：

- `-h`：同源基因名称文件
- `-n`：指定核酸序列文件
- `-a`：指定蛋白序列文件
- `-p`：线程数文件
- `-m`：指定比对工具（clusterw2、t_coffee、mafft、muscle 其中一个）
- `-g`：去除比对有 gap 的密码子
- `-k`：用 KaKs_Calculator 计算 KaKs 值
- `-o`：输出结果的目录
- `-f`：输出比对文件的格式

### 6.4 统计分析

上一步骤中计算 KaKs 输出目录中，包含两类文件，一个是包含比对计算的两个 Cdomain 的核酸序列、另一个是这两个 Cdomain 计算结果，其内容如下：

```tsv
Sequence	Method	Ka	Ks	Ka/Ks	P-Value(Fisher)	Length	S-Sites	N-Sites	Fold-Sites(0:2:4)	Substitutions	Syn-Subs	Nonsyn-Subs	Fold-Syn-Subs(0:2:4)	Fold-Nonsyn-Subs(0:2:4)	Divergence-Distance	Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)	GC(1:2:3)	ML-Score	AICc	Akaike-Weight	Model
Paenib_barc_KACC11450_GCF_013347305_1-Condensation_LCL-2371-2655-3.8-none-Paenib_dendri_2022CK_00834_GCF_029952845_1-Condensation_LCL-2400-2684-2.8-none	MA	0.0661559	0.337669	0.195919	1.9645e-17	852	224.486	627.514	NA	95	61.383	33.617	NA	NA	0.137695	4.04176:3.64164:2.62172:2.61015:1.52683:1	0.535211(0.637324:0.339789:0.628521)	-1481.94	NA	NA	NA
```

一共有两行，第一行是表头，第二行是具体数据，其中第一、四、五列分别为名称、Ks 值、Ka/Ks 比值。在这里，统计分析 Ks 值

#### 6.4.1 提取名称和 Ks 值

首先将结果文件中名称和 Ks 值提取出来生成新文件

```shell
cd ~/chenxy/Pae_rerun/result/KaKs
mkdir -p Ks_Cdomain

for location in $(cat 10_Location.lst); do
    ks_all="${location}";
    for i in $(ls ${location}_result_dir | grep ".kaks"); do
        ks=$(cat ${location}_result_dir/${i} | tail -n 1 | cut -f 4);
        ks_all="${ks_all},${ks}";
    done
    echo "${ks_all}" | sed "s/,/\n/g" >> Ks_Cdomain/${location}_ks.tsv
done

# 将十个文件的 Ks值 合并为一个文件
cd Ks_Cdomain
paste *.tsv > Ks_Cdomain_all.tsv

head -n 5 Ks_Cdomain_all.tsv
# 文件格式
# pmxA|DCL|Location1	pmxA|LCL|Location2	pmxA|LCL|Location3	pmxA|LCL|Location4	pmxE|C_Starter	pmxE|C|Location4	pmxE|LCL|Location2	pmxE|LCL|Location3	pmxE|LCL|Location5	pmxE|LCL|Location6
# 0.331552	0.337669	0.455732	0.403827	0.365208	0.359803	0.514165	0.45509	0.497861	0.428051
# 0.34714	0.358534	0.454202	0.470959	0.436422	0.34094	0.512858	0.473704	0.519772	0.408399
# 0.486849	0.529026	0.53781	0.429783	0.472441	0.484432	0.659364	0.505616	0.641966	0.507292
# 0.366262	0.380689	0.427723	0.367751	0.410639	0.34004	0.577856	0.407604	0.496559	0.45796

```

#### 6.4.2 统计频率

- R包统计

```R
setwd("~/chenxy/Pae_rerun/result/KaKs")

# 加载R包
library(pheatmap)
library(readr)
library(dplyr)

# define function: 将 Ks值[0,1.5]划分为不同区间，统计每个区间的Ks值的频率
bins_ks_freq <- function(column_data, breaks) {
    # 过滤NA
    column_data <- column_data[!is.na(column_data)]
    # 统计总个数
    total_count <- length(column_data)
    # 依据 Ks值的大小分组
    bins <- cut(column_data, breaks = breaks, include.lowest = TRUE)
    # 每个分组中元素的频数
    counts <- table(bins)
    counts_df <- as.data.frame(counts)
    colnames(counts_df) <- c("ks_bin", "count")
    # 每个分组中元素的频率 [0,100]
    counts_df$count_normolized <- (counts_df$count / total_count)*100
    return(counts_df)
}

# data
# in_df <- read_tsv("3search/6result/5ks/3paml/adomain_site.ks.tsv")
# in_df <- read_tsv("3search/6result/5ks/3paml/adomain_substrate.ks.tsv")
# in_df <- read_tsv("3search/6result/5ks/3paml/adomain_branch.ks.tsv")
# in_df <- read_tsv("3search/6result/5ks/3paml/tdomain_function.ks.tsv")
# in_df <- read_tsv("3search/6result/5ks/3paml/tdomain_site.ks.tsv")
# in_df <- read_tsv("3search/6result/5ks/3paml/cdomain_function.ks.tsv")
# in_df <- read_tsv("3search/6result/5ks/3paml/cdomain_site.ks.tsv")
# in_df <- read_tsv("3search/6result/5ks/3paml/pmx_site.ks.tsv")
in_df <- read_tsv("Ks_Cdomain/Ks_Cdomain_all.tsv")
## 将n/c值全部转变为NA
in_df[] <- lapply(in_df, function(x) as.numeric(as.character(x)))
## 删除全为NA的行
in_df <- in_df[rowSums(is.na(in_df)) < ncol(in_df), ]

# 打印 Ks值 的最小值、最大值和 >1.5的值
print(max(in_df, na.rm = TRUE))    # 1.79848
print(min(in_df, na.rm = TRUE))    #  0.00796442
print(in_df[in_df > 1.5], na.rm = TRUE)

# >1.5的值定义为1.5, 处理前过滤 NA
in_df[] <- lapply(in_df, function(x) {
  x[!is.na(x) & x > 1.5] <- 1.5
  return(x)
})

# 检查ks值分布 —— 是否存在小于0或大于1.5的值
all_in_range <- all(in_df >= 0 & in_df <= 1.5, na.rm = TRUE)
if (all_in_range) {
  print("YES")
} else {
  print("NO")
}

# breaks
breaks <- seq(0,1.5, length.out = 61)
ks_freq_list <- lapply(in_df, bins_ks_freq, breaks = breaks)

# combine
ks_bins <- ks_freq_list[[1]]$ks_bin
ks_freq_df <- data.frame(ks_bins = ks_bins)
for (i in seq_along(ks_freq_list)) {
    ks_freq_df[[names(in_df)[i]]] <- ks_freq_list[[i]]$count_normolized
}
ks_matrix <- as.matrix(ks_freq_df[, -1])
rownames(ks_matrix) <- ks_freq_df$ks_bin

# plot
plot <- pheatmap(
    t(ks_matrix),
    color = colorRampPalette(c("#ffffff", "#e62d2d"))(60),
    scale = "none",  # 不进行标准化处理
    cluster_rows = FALSE, cluster_cols = FALSE,
    border_color = NA,  # 设置单元格边框颜色
    fontsize_row = 7, fontsize_col = 7,
    main = "Ks Heatmap",
    # display_numbers = TRUE, fontsize_number = 5 # 显示单元格数值
    cellwidth = 12, cellheight = 18,  # 设置单元格大小
)

# "navy"
```
结果如图：
![freq_heatmap_Cdomain](/pic/freq_heatmap_Cdomain.png 'freq_heatmap_Cdomain')

- `nwr plot hh`统计

[nwr](https://github.com/wang-q/nwr)更新了绘制热图功能，在这里，可以使用另一种方法进行统计

`nwr plot hh`输入文件格式要求为两列，在第二列中分组
```tsv
Value	Group
0.127355	A1
0.136287	A1
0.151956	A1
0.123647	A2
0.126801	A2
0.153377	A2
```


```shell
cd ~/chenxy/Pae_rerun/result/KaKs

# 先生成 nwr plot hh需要的格式
for location in $(cat 10_Location.lst); do
    header=$(head -n 1 Ks_Cdomain/${location}_ks.tsv );
    tail -n +2 Ks_Cdomain/${location}_ks.tsv | 
        sed "s/$/\t$header/" >> Ks_Cdomain/${location}_ks.txt
done
cat Ks_Cdomain/*.txt >> Ks_Cdomain/Ks_Cdomain.tsv

# 生成热图
nwr plot hh Ks_Cdomain/Ks_Cdomain.tsv -g 2 --bins 40 --xl "" --yl "" --xmm 0,2  -o Ks_plot/Ks_Cdomain.tex 

# 由于名称中含下划线，打开tex文件将 yticklabels 中的 pmxE|C_Starter 修改为 pmxE|C\_Starter
# tex文件转化为pdf
tectonic Ks_plot/Ks_Cdomain.tex
magick -density 600 Ks_plot/Ks_Cdomain.pdf Ks_plot/Ks_Cdomain.png
```
参数：
> `-g`输入tsv的列数，一般为2
>
> `--bins`将横坐标区间拆分为多少小区间，默认40
>
> `--xl`/`--yl`横纵坐标的名称
>
> `--xmm`横坐标的最小值和最大值
>
> `-o`输出文件的名称，默认是stdout

结果如图：
![heatmap_Cdomain](/pic/Ks_Cdomain.png 'heatmap_Cdomain')