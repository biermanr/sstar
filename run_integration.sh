#!/bin/bash
set -euo pipefail

micromamba run -n sstar-build coverage run -m pytest -m integration
micromamba run -n sstar-build coverage erase
