# 本地下载 SRA 数据并上传服务器教程

服务器网络限制了对 NCBI/ENA 国际数据库的访问，需在本地（有网络的机器）下载后传至服务器。

## 需要下载的数据集

| 项目 | 样本数 | 说明 |
|------|--------|------|
| PRJNA542018 | ~40 | 16S V4，单端（merged PE） |
| PRJNA598825 | 30 | 16S，OLP 口腔微生物组 |
| PRJNA1049117 | 43 | 16S，OLP vs 健康对照 |
| PRJDB12280 | 77 | 16S，仅 Presort fraction（已筛选） |

服务器目标目录：`/home/usr/yel/exp-OLP_meta/pipeline/data/OLP/raw/meta/`

---

## 第一步：本地安装 iSeq

```bash
# 方式一：conda（推荐）
conda create -n iseq -c bioconda -c conda-forge iseq -y
conda activate iseq

# 方式二：直接下载脚本
wget https://raw.githubusercontent.com/BioOmics/iSeq/main/iseq
chmod +x iseq
```

验证安装：
```bash
iseq --version
```

---

## 第二步：下载各数据集

在本地建立工作目录：
```bash
mkdir -p ~/olp_download
cd ~/olp_download
```

### PRJNA542018（全量下载）

```bash
iseq -i PRJNA542018 -g -p 5 -o PRJNA542018/
```

### PRJNA598825（全量下载）

```bash
iseq -i PRJNA598825 -g -p 5 -o PRJNA598825/
```

### PRJNA1049117（按 accession list 下载）

项目号在 ENA/SRA 索引中缺失，需按 Run ID 列表下载。

新建文件 `PRJNA1049117_accessions.txt`，内容如下：
```
SRR27100551
SRR27100552
SRR27100553
SRR27100554
SRR27100555
SRR27100556
SRR27100557
SRR27100558
SRR27100559
SRR27100560
SRR27100561
SRR27100562
SRR27100563
SRR27100564
SRR27100565
SRR27100566
SRR27100567
SRR27100568
SRR27100569
SRR27100570
SRR27100571
SRR27100572
SRR27100573
SRR27100574
SRR27100575
SRR27100576
SRR27100577
SRR27100578
SRR27100579
SRR27100580
SRR27100581
SRR27100582
SRR27100583
SRR27100584
SRR27100585
SRR27100586
SRR27100587
SRR27100588
SRR27100589
SRR27100590
SRR27100591
SRR27100592
SRR27100593
```

下载：
```bash
mkdir -p PRJNA1049117
iseq -i PRJNA1049117_accessions.txt -g -p 5 -o PRJNA1049117/
```

### PRJDB12280（Presort fraction，筛选后 77 个样本）

从服务器取 accession list（服务器上已生成）：
```bash
scp yel@<SERVER_IP>:/home/usr/yel/exp-OLP_meta/pipeline/data/OLP/raw/meta/PRJDB12280/accession_list.txt \
    ./PRJDB12280_accessions.txt
```

下载：
```bash
mkdir -p PRJDB12280
iseq -i PRJDB12280_accessions.txt -g -p 5 -o PRJDB12280/
```

---

## 第三步：验证下载完整性

```bash
# 列出各项目文件数量
for proj in PRJNA542018 PRJNA598825 PRJNA1049117 PRJDB12280; do
    count=$(ls ${proj}/*.fastq.gz 2>/dev/null | wc -l)
    echo "${proj}: ${count} 个 .fastq.gz 文件"
done
```

如有未下载完的，重新运行同一命令即可（iSeq 自动跳过已存在文件）。

---

## 第四步：上传到服务器

### 方式A：scp（简单，适合文件较少）

```bash
SERVER="yel@<SERVER_IP>"
REMOTE="/home/usr/yel/exp-OLP_meta/pipeline/data/OLP/raw/meta"

scp -r PRJNA542018/  ${SERVER}:${REMOTE}/PRJNA542018/
scp -r PRJNA598825/  ${SERVER}:${REMOTE}/PRJNA598825/
scp -r PRJNA1049117/ ${SERVER}:${REMOTE}/PRJNA1049117/
scp -r PRJDB12280/   ${SERVER}:${REMOTE}/PRJDB12280/
```

### 方式B：rsync（推荐，支持断点续传）

```bash
SERVER="yel@<SERVER_IP>"
REMOTE="/home/usr/yel/exp-OLP_meta/pipeline/data/OLP/raw/meta"

for proj in PRJNA542018 PRJNA598825 PRJNA1049117 PRJDB12280; do
    echo "上传 ${proj}..."
    rsync -avz --progress ${proj}/ ${SERVER}:${REMOTE}/${proj}/
done
```

> **注意**：将 `<SERVER_IP>` 替换为实际服务器 IP 或主机名。

---

## 第五步：验证服务器上的文件

上传完成后在服务器上确认：

```bash
for proj in PRJNA542018 PRJNA598825 PRJNA1049117 PRJDB12280; do
    count=$(ls /home/usr/yel/exp-OLP_meta/pipeline/data/OLP/raw/meta/${proj}/*.fastq.gz 2>/dev/null | wc -l)
    echo "${proj}: ${count} 个文件"
done
```

---

## 下一步

数据上传完成后，在服务器上运行 QIIME2 处理：

```bash
conda activate qiime2-amplicon-2023
cd /home/usr/yel/exp-OLP_meta/pipeline

# PRJNA542018：单端 V4，需 cutadapt 截除引物
nohup bash scripts/11_qiime2_PRJNA542018.sh > logs/qiime2_PRJNA542018.log 2>&1 &

# CRA008410：双端 V3-V4
nohup bash scripts/12_qiime2_CRA008410.sh > logs/qiime2_CRA008410.log 2>&1 &

# PRJNA598825、PRJNA1049117、PRJDB12280：QIIME2 脚本待新建（需确认引物和测序类型）
```
