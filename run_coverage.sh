#!/bin/bash
set -euo pipefail

micromamba run -n sstar-build coverage run -m pytest -m "not integration"
micromamba run -n sstar-build coverage combine
micromamba run -n sstar-build coverage report -m
micromamba run -n sstar-build coverage erase
