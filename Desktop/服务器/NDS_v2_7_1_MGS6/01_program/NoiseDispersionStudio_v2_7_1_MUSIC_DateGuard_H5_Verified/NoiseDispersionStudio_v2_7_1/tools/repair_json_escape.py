from __future__ import annotations
from pathlib import Path
import argparse
import json
import re
import shutil


def repair_json_text(text: str) -> str:
    token_re = re.compile(r'"(?:\\.|[^"\\])*"')

    def fix_token(match: re.Match[str]) -> str:
        token = match.group(0)
        body = token[1:-1]
        is_drive_path = bool(re.match(r'^[A-Za-z]:\\', body))
        is_unc_path = body.startswith('\\\\')
        if is_drive_path or is_unc_path:
            if is_unc_path and not body.startswith(r'\\\\'):
                body=r'\\'+body
            out = []
            i = 0
            while i < len(body):
                ch = body[i]
                if ch != '\\':
                    out.append(ch)
                    i += 1
                    continue
                if i + 1 < len(body) and body[i + 1] == '\\':
                    out.append('\\\\')
                    i += 2
                    continue
                if i + 1 < len(body) and body[i + 1] == '"':
                    out.append('\\"')
                    i += 2
                    continue
                out.append('\\\\')
                i += 1
            return '"' + ''.join(out) + '"'
        return token

    return token_re.sub(fix_token, text)


def repair_file(path: Path) -> str:
    raw = path.read_text(encoding='utf-8-sig')
    try:
        json.loads(raw)
        return 'OK'
    except json.JSONDecodeError as first:
        repaired = repair_json_text(raw)
        try:
            data = json.loads(repaired)
        except json.JSONDecodeError:
            return f'SKIP: not a Windows-path escape problem ({first})'
        backup = path.with_name(path.name + '.invalid_escape_backup')
        if not backup.exists():
            shutil.copy2(path, backup)
        temp = path.with_name(path.name + '.repairing')
        temp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
        temp.replace(path)
        return f'REPAIRED (backup: {backup.name})'


def candidates(root: Path):
    direct = [
        root / 'hdf_import_audit.json',
        root / 'project.json',
        root / 'config.json',
        root / 'mapping.json',
        root / 'phase_settings.json',
        root / 'selection.json',
        root / 'last_run.json',
    ]
    out = root / 'output'
    if out.exists():
        direct += sorted(out.glob('*.json'))
        phase = out / 'phase_guided'
        if phase.exists():
            direct += sorted(phase.glob('*.json'))
    seen = set()
    for p in direct:
        if p.exists() and p not in seen:
            seen.add(p)
            yield p


def main() -> int:
    ap = argparse.ArgumentParser(description='Repair legacy invalid Windows backslashes in Studio JSON files.')
    ap.add_argument('project', type=Path, help='Studio project directory (contains project.json)')
    a = ap.parse_args()
    root = a.project.resolve()
    if not (root / 'project.json').exists():
        raise SystemExit(f'Not a Studio project directory: {root}')
    repaired = 0
    for p in candidates(root):
        result = repair_file(p)
        print(f'{p}: {result}')
        repaired += result.startswith('REPAIRED')
    print(f'Finished. Repaired {repaired} file(s).')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
