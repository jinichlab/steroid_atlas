#!/bin/bash
# Launch the Steroid Atlas visualizer.
# Adjust PORT and LD_LIBRARY_PATH if needed for your environment.

set -e
cd "$(dirname "$0")/.."

PORT="${PORT:-2730}"
HOST="${HOST:-0.0.0.0}"

# If marimo/rdkit can't load system libstdc++ (common on servers with out-of-date
# system libs), point at a newer libstdc++ from a conda env. Comment out or
# override if not needed.
if [ -d "$HOME/miniconda3/lib" ]; then
  export LD_LIBRARY_PATH="$HOME/miniconda3/lib:${LD_LIBRARY_PATH:-}"
fi

# Increase marimo's per-output byte cap so large protein tables don't get truncated
export MARIMO_OUTPUT_MAX_BYTES="${MARIMO_OUTPUT_MAX_BYTES:-50000000}"

echo "Launching Steroid Atlas visualizer on ${HOST}:${PORT}"
echo "From your laptop, open an SSH tunnel first:"
echo "  ssh -L ${PORT}:localhost:${PORT} -o ServerAliveInterval=30 <you>@<host>"
echo "Then browse to http://localhost:${PORT}"
echo ""

exec marimo run app/visualizer.py --port "$PORT" --host "$HOST"
