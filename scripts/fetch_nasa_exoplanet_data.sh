#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
RAW_DIR="$ROOT_DIR/data/raw"

mkdir -p "$RAW_DIR"

PS_URL="https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,hostname,disc_year,discoverymethod,pl_orbper,pl_rade,pl_insol,pl_eqt,pl_trandep,pl_trandur,st_teff,st_rad,st_mass,st_met,sy_dist,sy_jmag,sy_hmag,sy_kmag+from+pscomppars+where+tran_flag=1+and+pl_rade+between+0.5+and+6+and+pl_orbper+between+0.3+and+100&format=csv"
TOI_URL="https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+toi,toipfx,tfopwg_disp,pl_pnum,st_tmag,ra,dec,pl_orbper,pl_trandurh,pl_trandep,pl_rade,pl_insol,pl_eqt,st_teff,st_rad+from+toi+where+tfopwg_disp=%27PC%27+and+pl_rade+between+0.5+and+6+and+pl_orbper+between+0.3+and+100&format=csv"

echo "Fetching PSCompPars snapshot..."
curl -L "$PS_URL" -o "$RAW_DIR/pscomppars_small_transit.csv"

echo "Fetching TOI snapshot..."
curl -L "$TOI_URL" -o "$RAW_DIR/toi_pc_small.csv"

echo "Finished. Files written to $RAW_DIR"
