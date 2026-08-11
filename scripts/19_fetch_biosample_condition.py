#!/usr/bin/env python3
"""
BioSample 属性查询脚本
对 combined_metadata.tsv 中 condition=FILL_ME 的样本，
批量查询 NCBI BioSample API，自动推断 OLP/Control 分组。

用法：
    python3 scripts/19_fetch_biosample_condition.py \
        --input  results/OLP/metadata/combined_metadata.tsv \
        --output results/OLP/metadata/combined_metadata_annotated.tsv \
        --log    results/OLP/metadata/biosample_fetch.log

选项：
    --dry-run   只打印查询结果，不写入文件
    --delay     每次请求间隔秒数（默认 0.34，约 3 req/s，NCBI 限速）
    --email     提供给 NCBI 的联系邮箱（建议填写）
"""

import argparse
import csv
import sys
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict

# ── OLP / Control 关键词 ─────────────────────────────────────────────────────
OLP_KEYWORDS  = [
    "olp", "oral lichen planus", "lichen planus", "lp patient",
    "patient", "case", "diseased", "affected",
]
CTRL_KEYWORDS = [
    "control", "ctrl", "healthy", "normal", "hc", "non-olp",
    "unaffected", "without", "disease free",
]


def infer_from_text(text: str) -> str:
    if not text:
        return ""
    t = text.lower()
    if any(k in t for k in OLP_KEYWORDS):
        return "OLP"
    if any(k in t for k in CTRL_KEYWORDS):
        return "Control"
    return ""


# ── NCBI eUtils ──────────────────────────────────────────────────────────────
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

# BioSample 属性名（优先级顺序）
CONDITION_ATTRS = [
    "disease_state", "disease", "health_state", "host_health_state",
    "condition", "status", "diagnosis", "phenotype",
    "host_disease", "isolation_source", "source_name",
    "patient_group", "group", "sample_type", "tissue",
    "subject_is_affected",
]


def fetch_biosample_xml(biosample_id: str, email: str, delay: float) -> ET.Element | None:
    """查询单个 BioSample，返回 XML Element 或 None"""
    params = {
        "db": "biosample",
        "id": biosample_id,
        "retmode": "xml",
    }
    if email:
        params["email"] = email

    url = EFETCH_URL + "?" + urllib.parse.urlencode(params)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "OLP-metaanalysis/1.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = resp.read()
        time.sleep(delay)
        return ET.fromstring(data)
    except Exception as e:
        return None


def parse_biosample_attrs(root: ET.Element) -> dict:
    """从 BioSample XML 提取所有属性"""
    attrs = {}
    # BioSample 属性在 <Attributes> 下的 <Attribute>
    for attr in root.findall(".//{http://www.ncbi.nlm.nih.gov/biosample}Attribute"):
        name = (attr.get("attribute_name") or attr.get("harmonized_name") or "").lower().strip()
        val  = (attr.text or "").strip()
        if name and val:
            attrs[name] = val
    # 兼容无命名空间的格式
    for attr in root.findall(".//Attribute"):
        name = (attr.get("attribute_name") or attr.get("harmonized_name") or "").lower().strip()
        val  = (attr.text or "").strip()
        if name and val:
            attrs.setdefault(name, val)
    # 样本标题
    for tag in ["{http://www.ncbi.nlm.nih.gov/biosample}Description",
                "{http://www.ncbi.nlm.nih.gov/biosample}Title", "Description", "Title"]:
        el = root.find(f".//{tag}")
        if el is not None and el.text:
            attrs.setdefault("_title", el.text.strip())
    return attrs


def infer_condition_from_attrs(attrs: dict) -> tuple[str, str]:
    """
    返回 (condition, evidence)
    condition: OLP / Control / ""（无法推断）
    evidence:  依据描述
    """
    # 1. 明确属性优先
    for key in CONDITION_ATTRS:
        val = attrs.get(key, "")
        if not val:
            continue
        cond = infer_from_text(val)
        if cond:
            return cond, f"{key}={val!r}"

    # 2. 所有属性值全文搜索
    for key, val in attrs.items():
        if key.startswith("_"):
            continue
        cond = infer_from_text(val)
        if cond:
            return cond, f"{key}={val!r}"

    # 3. 标题
    title_cond = infer_from_text(attrs.get("_title", ""))
    if title_cond:
        return title_cond, f"title={attrs['_title']!r}"

    return "", "no_match"


# ── 主流程 ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",   default="results/OLP/metadata/combined_metadata.tsv")
    parser.add_argument("--output",  default="results/OLP/metadata/combined_metadata_annotated.tsv")
    parser.add_argument("--log",     default="results/OLP/metadata/biosample_fetch.log")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--delay",   type=float, default=0.34,
                        help="API 请求间隔秒数（默认 0.34 ≈ 3 req/s）")
    parser.add_argument("--email",   default="",
                        help="提供给 NCBI 的联系邮箱（建议填）")
    parser.add_argument("--study",   default="",
                        help="只处理指定数据集（如 PRJNA542018），留空处理所有 FILL_ME")
    args = parser.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        sys.exit(f"[ERROR] 找不到输入文件：{in_path}")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    # 读取现有 metadata
    with open(in_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = list(reader)

    # 找需要查询的行
    to_fetch = [
        r for r in rows
        if r.get("condition", "").strip() == "FILL_ME"
        and r.get("biosample_id", "").strip()
        and (not args.study or r.get("study_id") == args.study)
    ]

    # 按 biosample_id 去重（同一 biosample 可能有多个 run）
    biosample_ids = list({r["biosample_id"].strip() for r in to_fetch})

    print(f"[INFO] 需要查询 BioSample: {len(biosample_ids)} 个（覆盖 {len(to_fetch)} 条记录）")
    if not biosample_ids:
        print("[OK] 没有需要补全的样本，退出。")
        return

    # 查询 BioSample
    biosample_cache: dict[str, dict] = {}   # biosample_id -> attrs
    log_lines = []

    for i, bio_id in enumerate(biosample_ids, 1):
        print(f"  [{i:3d}/{len(biosample_ids)}] {bio_id} ...", end=" ", flush=True)
        root = fetch_biosample_xml(bio_id, args.email, args.delay)
        if root is None:
            print("FAILED")
            log_lines.append(f"{bio_id}\tFAILED\t\t")
            biosample_cache[bio_id] = {}
            continue

        attrs = parse_biosample_attrs(root)
        cond, evidence = infer_condition_from_attrs(attrs)
        biosample_cache[bio_id] = {"cond": cond, "evidence": evidence, "attrs": attrs}

        status = cond if cond else "UNKNOWN"
        print(f"{status:8s}  ({evidence})")
        log_lines.append(f"{bio_id}\t{status}\t{evidence}\t"
                         + "; ".join(f"{k}={v}" for k, v in list(attrs.items())[:8]
                                     if not k.startswith("_")))

    # 统计
    n_olp  = sum(1 for v in biosample_cache.values() if v.get("cond") == "OLP")
    n_ctrl = sum(1 for v in biosample_cache.values() if v.get("cond") == "Control")
    n_unk  = sum(1 for v in biosample_cache.values() if not v.get("cond"))
    print(f"\n[结果] OLP={n_olp}  Control={n_ctrl}  UNKNOWN={n_unk}")

    if args.dry_run:
        print("[dry-run] 不写入文件。")
        return

    # 更新 condition 字段
    updated = 0
    still_fill = 0
    for row in rows:
        if row.get("condition", "").strip() != "FILL_ME":
            continue
        bio_id = row.get("biosample_id", "").strip()
        if not bio_id or bio_id not in biosample_cache:
            still_fill += 1
            continue
        cond = biosample_cache[bio_id].get("cond", "")
        if cond:
            row["condition"] = cond
            updated += 1
        else:
            still_fill += 1
            # 把找到的属性值写入 notes，方便手工判断
            attrs = biosample_cache[bio_id].get("attrs", {})
            hints = "; ".join(
                f"{k}={v}" for k, v in attrs.items()
                if k in CONDITION_ATTRS and v and not k.startswith("_")
            )
            if hints:
                row["notes"] = hints

    # 写输出
    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # 写日志
    with open(args.log, "w", encoding="utf-8") as f:
        f.write("biosample_id\tcondition\tevidence\tattributes_preview\n")
        f.write("\n".join(log_lines))

    print(f"\n[OK] 已更新 {updated} 条记录 → {args.output}")
    if still_fill > 0:
        print(f"[!]  仍有 {still_fill} 条无法自动判断，condition 保持 FILL_ME")
        print(f"     请查看 notes 列的属性提示，手工编辑 {args.output}")
    print(f"[OK] 查询日志 → {args.log}")


if __name__ == "__main__":
    main()
