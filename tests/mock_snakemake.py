"""
模拟 Snakemake 注入的 snakemake 对象，用于独立测试各模块脚本。
用法：在测试脚本中 import 并设置全局变量后，再 exec() 目标脚本。
"""

class _AttrDict:
    """
    支持三种访问方式：
      obj.key          → 具名访问
      obj["key"]       → 字符串键访问
      obj[0]           → 整数位置索引（Snakemake 的 input[0]/output[0] 风格）
    """
    def __init__(self, d=None):
        self._d = d or {}

    def __getattr__(self, name):
        if name.startswith("_"):
            return super().__getattribute__(name)
        return self._d[name]

    def __getitem__(self, key):
        if isinstance(key, int):
            # 整数索引：按插入顺序取第 key 个值
            return list(self._d.values())[key]
        return self._d[key]

    def __iter__(self):
        return iter(self._d.values())

    def __repr__(self):
        return f"AttrDict({self._d})"


class MockSnakemake:
    def __init__(self, input=None, output=None, params=None, wildcards=None):
        self.input     = _AttrDict(input  or {})
        self.output    = _AttrDict(output or {})
        self.params    = _AttrDict(params or {})
        self.wildcards = _AttrDict(wildcards or {})


def make_snakemake(input=None, output=None, params=None, wildcards=None):
    return MockSnakemake(input=input, output=output, params=params, wildcards=wildcards)
