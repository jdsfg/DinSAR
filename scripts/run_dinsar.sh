#!/bin/bash
# 兼容入口: 实际流程只维护仓库根目录的 run_dinsar.sh, 这里直接转交。
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/run_dinsar.sh" "$@"
