#!/bin/sh
# Upload the web app to the Hugging Face Space, which rebuilds it from the
# Dockerfile. Needs `hf auth login` with write access to the Space.
#
#   sh scripts/deploy_space.sh                      # jinichlab/steroid-atlas
#   SPACE=you/your-space sh scripts/deploy_space.sh
#
# Only what the Docker build needs is uploaded — never .env*, node_modules or
# the RAG index (the Space downloads that from its private dataset at startup;
# see scripts/start_space.sh). The analysis side (data/, analysis/, literature/,
# marimo/, …) stays in git only.
#
# The Space is made an exact mirror of what's uploaded: `--delete "*"` removes
# any remote file not in this upload (Hugging Face never deletes .gitattributes),
# so files from old layouts or a mistaken upload can't linger.
set -eu
cd "$(dirname "$0")/.."
SPACE="${SPACE:-jinichlab/steroid-atlas}"

# Refuse to run from anything but the steroid_atlas repo root.
for f in Dockerfile web/package.json rag/rag_server.py scripts/start_space.sh space/README.md; do
  [ -e "$f" ] || { echo "deploy_space.sh: $f not found — not the steroid_atlas repo root" >&2; exit 1; }
done

# Stage exactly what the Docker build needs. The repo README is for GitHub;
# the Space gets space/README.md, whose front matter configures the Space.
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
rsync -a \
  --exclude node_modules --exclude .next --exclude '.env*' --exclude '*.tsbuildinfo' \
  --exclude __pycache__ \
  web rag "$STAGE"/
mkdir -p "$STAGE/scripts"
cp Dockerfile .dockerignore LICENSE "$STAGE"/
cp scripts/start_space.sh "$STAGE/scripts/"
cp space/README.md "$STAGE/README.md"

hf upload "$SPACE" "$STAGE" . --repo-type space \
  --delete "*" \
  --commit-message "${MSG:-Deploy from $(git rev-parse --short HEAD 2>/dev/null || echo local)}"
