#!/usr/bin/env bash
set -e
cd -- "$(dirname -- "$0")"
if [[ -x .venv/bin/python ]]; then
  .venv/bin/python run.py download
else
  python3 run.py download
fi
