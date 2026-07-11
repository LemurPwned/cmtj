#!/usr/bin/env python3
"""Claude Code PostToolUse hook wrapper around check_units.py.

Reads the standard hook JSON payload from stdin, and if the edited/written
file is a .py containing "cmtj", runs check_units.py on it and prints any
warnings. Always exits 0 -- informational only, never blocks the tool call.

Wire it up in .claude/settings.json (project or user level):

{
  "hooks": {
    "PostToolUse": [
      {
        "matcher": "Write|Edit",
        "hooks": [
          {
            "type": "command",
            "command": "python3 skills/cmtj-simulation/scripts/post_edit_hook.py"
          }
        ]
      }
    ]
  }
}
"""

import json
import os
import sys

from check_units import check


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except json.JSONDecodeError:
        return 0

    file_path = payload.get("tool_input", {}).get("file_path", "")
    if not file_path.endswith(".py") or not os.path.exists(file_path):
        return 0

    with open(file_path) as f:
        if "cmtj" not in f.read():
            return 0

    check(file_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
