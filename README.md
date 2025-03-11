# 同义/非同义突变计算
- 定义：
> 同义突变率ds/Ks：发生同义突变的位点数目（Ds）与可能发生同义突变的所有位点总数的比值
> 
> 非同义突变率dN/Ka：发生非同义突变的位点数目（Dn）与可能发生非同义突变的所有位点总数的比值

- 算法模型：
> 近似法：
> 最大似然法：利用概率论完成近似法的三步骤。代表模型YN。
> 计数法：
>

## 0 软件
### antismash  v7.1.0
- 安装`conda`
```shell
# 下载安装conda
mkdir ~/miniconda3
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
- 使用`conda`手动安装antismash
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

# 用法
conda activate antismash
antismash my_input.gbk
```

### KaKs_Calculator v3.0
[下载地址](https://github.com/Chenglin20170390/KaKs_Calculator-3.0.git)
```shell
cd ~/biosoft
git clone https://github.com/Chenglin20170390/KaKs_Calculator-3.0.git
cd KaKs_Calculator-3.0/bin
make
echo 'export PATH=$HOME/biosoft/KaKs_Calculator-3.0/bin:$PATH'  >> $HOME/.bashrc
source ~/.bashrc
```

### ParaAT  v2.0
整合了计算KaKs所需的一整套分析，包括蛋白序列比对、根据比对结果回译成codon对应的核酸比对结果、计算KaKs值（使用KaKs_Calculator）   
[下载地址1](https://ngdc.cncb.ac.cn/tools/paraat)   
[下载地址2](https://github.com/jdebarry/paraat.git)   
```shell
# 解压
unzip ParaAT2.0.tar.gz
echo 'export PATH=$HOME/biosoft/ParaAT2.0:$PATH'  >> $HOME/.bashrc
source ~/.bashrc

# 依赖工具下载（蛋白比对）
brew install mafft

# 测试是否安装成功
ParaAT.pl -h test.homologs -n test.cds -a test.pep -p proc -m mafft -f axt -g -k -o result_dir
# 成功安装后有输出目录result_dir
```

#### ParaAT2.0用法
```shell
ParaAT.pl -h test.homologs -n test.cds -a test.pep -p proc -m mafft -f axt -g -k -o result_dir
```
- `-h`：同源基因名称文件
- `-n`：指定核酸序列文件
- `-a`：指定蛋白序列文件
- `-p`：线程数
- `-m`：指定比对工具（clusterw2、t_coffee、mafft、muscle其中一个）
- `-g`：去除比对有gap的密码子
- `-k`：用KaKs_Calculator计算KaKs值
- `-o`：输出结果的目录
- `-f`：输出比对文件的格式


## 1 数据下载





## 2 识别直系同源基因





## 3 直系同源基因比对





## 4 估算系统发育树




## 5 