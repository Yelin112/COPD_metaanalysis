"""
步骤4：解析 MetaCyc reactions.tab，提取酶-代谢反应信息
输出：microbial_metabolic_reactions.tsv
  列：ec | left_compounds | right_compounds | direction
对应原 Perl 脚本：4_extract_metacyc_reaction_info.pl

Snakemake 注入变量：
  snakemake.input[0]  — MetaCyc reactions.tab 路径
  snakemake.output[0] — 输出 TSV
"""

import re
import pandas as pd
from collections import defaultdict


def parse_reactions_tab(filepath: str) -> list[dict]:
    """
    MetaCyc reactions.tab 是多行块格式：
    每个反应块以 'UNIQUE-ID' 开头，以 '//' 结束
    相关字段：EC-NUMBER / LEFT / RIGHT / REACTION-DIRECTION
    """
    records = []
    current = defaultdict(list)

    with open(filepath, encoding="latin-1") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("//"):
                if current.get("EC-NUMBER") and (current.get("LEFT") or current.get("RIGHT")):
                    for ec in current["EC-NUMBER"]:
                        ec = ec.strip()
                        if not ec or ec == "?":
                            continue
                        records.append({
                            "ec":              ec,
                            "left_compounds":  ";".join(current.get("LEFT", [])),
                            "right_compounds": ";".join(current.get("RIGHT", [])),
                            "direction":       current.get("REACTION-DIRECTION", [""])[0].strip(),
                        })
                current = defaultdict(list)
                continue

            if " - " in line:
                key, _, value = line.partition(" - ")
                key = key.strip()
                value = value.strip()
                if key in ("EC-NUMBER", "LEFT", "RIGHT", "REACTION-DIRECTION"):
                    current[key].append(value)

    return records


def main():
    records = parse_reactions_tab(snakemake.input[0])
    df = pd.DataFrame(records)
    df.to_csv(snakemake.output[0], sep="\t", index=False)


main()
