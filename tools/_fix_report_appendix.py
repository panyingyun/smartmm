"""Rebuild research/03-ball-parameterization.md appendix cleanly.

The report body was authored via the write/edit tools (valid UTF-8).  Only the
appended appendix got mangled (a UTF-16 log was decoded as UTF-8).  This script
truncates at the last known-good body marker and re-appends a clean appendix.
"""
import os

RES = r"E:\panyingyun\smartmm\research"
REPORT = os.path.join(RES, "03-ball-parameterization.md")

# last line of the authored body
BODY_TAIL = "3. ❌ **不要**用 `K[idx][:,idx]` 去验证单元矩阵，也不要用 `tile` 写 scatter 索引。".encode("utf-8")

b = open(REPORT, "rb").read()
i = b.find(BODY_TAIL)
if i < 0:
    raise SystemExit("body tail marker not found -- inspect the file manually")
head = b[: i + len(BODY_TAIL)]
head.decode("utf-8")  # must not raise
print(f"body kept: {len(head)} bytes, valid UTF-8")

sweep = open(os.path.join(RES, "_sweep_convexity_log.txt"), encoding="utf-8").read()

import re
SKIP = re.compile(
    r"seconds:|Delauniz|Recover|Refin|Smooth|Improv|Writing|parsing switches|"
    r"^success$|Statistics:|^Input |^Mesh |Steiner|^Output |^Total "
)
lines = [l for l in sweep.splitlines() if l.strip() and not SKIP.search(l)]
clean = "\n".join(lines)

appendix = (
    "\n\n---\n\n"
    "## 附录 A：非凸度扫描的完整实测数据\n\n"
    "`tools/sweep_convexity.py` 的输出（已剔除 TetGen 命令行噪声）。\n"
    "列含义：`sphArea/4pi` = 边界球面映射的有符号立体角总和 / 4π（**1.0 = 双射**）；\n"
    "`A:*` = 普通 P1 harmonic 映射到球；`B:*` = 同一（径向投影）边界下的 VSEM；\n"
    "`B:gap` = $E_V(f) - \\frac32 V(f_{\\text{image}})$。\n\n"
    "```\n" + clean + "\n```\n"
)

out = head + appendix.encode("utf-8")
open(REPORT, "wb").write(out)
check = open(REPORT, "rb").read()
assert check.count(b"\x00") == 0, "NUL bytes still present"
check.decode("utf-8")
print(f"report written: {len(check)} bytes, 0 NUL bytes, valid UTF-8")
