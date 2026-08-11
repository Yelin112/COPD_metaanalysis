#!/usr/bin/env python3
"""
诊断脚本：打印 BioSample 的所有原始属性，
用于确认实际字段名，再调整 script 19 的关键词映射。

用法：
    python3 scripts/19b_debug_biosample.py SAMN11606556 SAMN11606540
    python3 scripts/19b_debug_biosample.py --from-tsv results/OLP/metadata/combined_metadata.tsv --study PRJNA542018 --n 3
"""
from __future__ import annotations

import argparse
import sys
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch_xml(biosample_id: str, email: str = "") -> str | None:
    params = {"db": "biosample", "id": biosample_id, "retmode": "xml"}
    if email:
        params["email"] = email
    url = EFETCH_URL + "?" + urllib.parse.urlencode(params)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "OLP-debug/1.0"})
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.read().decode("utf-8")
    except Exception as e:
        return None


def dump_attrs(xml_str: str, biosample_id: str):
    print(f"\n{'='*60}")
    print(f"BioSample: {biosample_id}")
    print(f"{'='*60}")
    try:
        root = ET.fromstring(xml_str)
    except ET.ParseError as e:
        print(f"  [XML解析失败] {e}")
        print(xml_str[:500])
        return

    # 标题 / 描述
    for tag in ["Title", "Description"]:
        for el in root.iter():
            if el.tag.split("}")[-1] == tag and el.text:
                print(f"  {tag}: {el.text.strip()}")

    # 所有 Attribute
    attrs_found = []
    for attr in root.iter():
        if attr.tag.split("}")[-1] == "Attribute":
            name = attr.get("attribute_name") or attr.get("harmonized_name") or ""
            val  = (attr.text or "").strip()
            attrs_found.append((name, val))

    if attrs_found:
        print(f"\n  Attributes ({len(attrs_found)} 个):")
        for name, val in attrs_found:
            print(f"    {name!r:40s} = {val!r}")
    else:
        print("\n  [未找到 Attribute 节点，打印原始 XML 前 800 字符]")
        print(xml_str[:800])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ids", nargs="*", help="BioSample ID（如 SAMN11606556）")
    parser.add_argument("--from-tsv", default="", help="从 combined_metadata.tsv 读取 ID")
    parser.add_argument("--study", default="", help="配合 --from-tsv 过滤数据集")
    parser.add_argument("--n", type=int, default=3, help="最多查几个（默认 3）")
    parser.add_argument("--email", default="")
    args = parser.parse_args()

    ids = list(args.ids)

    if args.from_tsv:
        import csv
        with open(args.from_tsv, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            seen = set()
            for row in reader:
                if args.study and row.get("study_id") != args.study:
                    continue
                if row.get("condition", "").strip() == "FILL_ME":
                    bid = row.get("biosample_id", "").strip()
                    if bid and bid not in seen:
                        seen.add(bid)
                        ids.append(bid)
                        if len(ids) >= args.n:
                            break

    if not ids:
        sys.exit("请提供 BioSample ID 或 --from-tsv")

    for i, bid in enumerate(ids[:args.n]):
        xml_str = fetch_xml(bid, args.email)
        if xml_str is None:
            print(f"\n[FAILED] 无法获取 {bid}")
        else:
            dump_attrs(xml_str, bid)
        if i < len(ids) - 1:
            time.sleep(0.4)


if __name__ == "__main__":
    main()
