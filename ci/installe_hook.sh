#!/usr/bin/env bash
# Branche la CI locale sur `git push` (versionnée dans ci/hooks, pas dans .git/hooks).
cd "$(dirname "$0")/.." && git config core.hooksPath ci/hooks && echo "pre-push branché : ci/local.sh juge chaque push"
