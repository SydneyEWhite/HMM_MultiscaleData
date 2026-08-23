#!/bin/bash

set -euo pipefail

mkdir -p logs sim_results

jid=$(sbatch --parsable submit_array.sh)
echo "array job:   $jid"

cid=$(sbatch --parsable --dependency=afterany:"$jid" submit_collect.sh)
echo "collect job: $cid"
