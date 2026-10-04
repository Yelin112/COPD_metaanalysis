#!/bin/bash
# HUMAnN3 (v3.9) 数据库下载脚本：镜像优先，wget -c 断点续传，下载后自动解压
# 用法: bash scripts/31_download_humann_db.sh <utility_mapping|chocophlan|uniref90>
#
# 版本依据: humann v3.9 源码 humann_databases.py 的 URL 表
#   - chocophlan full  = full_chocophlan.v201901_v31.tar.gz   (官方源，镜像无此版本)
#   - uniref90_diamond = uniref90_annotated_v201901b_full.tar.gz (unitn 镜像有)
#   - utility_mapping  = full_mapping_v201901b.tar.gz        (unitn 镜像有)
# 注意: 此脚本只负责下载解压；装好 humann 后再用 humann_config --update 指向目录。
set -u

DB=${1:-utility_mapping}
DEST=${DEST:-/home/exs/data3/Software_Database}
MIRROR="http://cmprod1.cibio.unitn.it/databases/HUMAnN"
HUTT="http://huttenhower.sph.harvard.edu/humann_data"
mkdir -p "$DEST"

case $DB in
  utility_mapping)
    URLS=("$MIRROR/full_mapping_v201901b.tar.gz"
          "$HUTT/full_mapping_v201901b.tar.gz")
    SUBDIR=utility_mapping ;;
  chocophlan)
    URLS=("$HUTT/chocophlan/full_chocophlan.v201901_v31.tar.gz")
    SUBDIR=chocophlan ;;
  uniref90)
    URLS=("$MIRROR/uniref90_annotated_v201901b_full.tar.gz"
          "$HUTT/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz")
    SUBDIR=uniref ;;
  *)
    echo "未知数据库: $DB (可选 utility_mapping|chocophlan|uniref90)"; exit 1 ;;
esac

# 逐个源尝试，wget -c 支持断点续传（中断后重跑本脚本即可）
OK=0
for url in "${URLS[@]}"; do
  fname=$(basename "$url")
  echo ">> 尝试: $url"
  if wget -c --tries=3 --timeout=60 -O "$DEST/$fname" "$url"; then OK=1; break; fi
done
[ "$OK" = 1 ] || { echo "所有源都下载失败"; exit 1; }

# 解压；若 tarball 只有一层顶层目录则上移拍平
mkdir -p "$DEST/$SUBDIR"
tar -xzf "$DEST/$fname" -C "$DEST/$SUBDIR" || { echo "解压失败"; exit 1; }
top=$(ls -A "$DEST/$SUBDIR")
if [ "$(echo "$top" | wc -l)" = 1 ] && [ -d "$DEST/$SUBDIR/$top" ]; then
  mv "$DEST/$SUBDIR/$top"/* "$DEST/$SUBDIR/" && rmdir "$DEST/$SUBDIR/$top"
fi
echo "完成: $DEST/$SUBDIR  (tarball 已保留在 $DEST/，确认无误后可删)"
