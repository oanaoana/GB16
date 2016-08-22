#!/home/maxhutch/anaconda3/bin/python3

import json
from os import getcwd

from nekpy.dask.subgraph import series
from nekpy.dask.tasks import configure, prepare, run
from nekpy.dask.utils import outer_product, work_name
from nekpy.config import config as cfg

from nekpy.dask import run_all

from sys import argv

nekmpi_path = cfg.nekmpi

with open(argv[1], "r") as f:
    base = json.load(f)

with open(argv[2], "r") as f:
    tusr = f.read()

base["prefix"] = 'swp' 

# Take simple outer product of contents of sweep file
candidates = [{"nodes" : 2**x} for x in range(3)]

# Filter out the cases we don't want
overrides = []
for c in candidates:
    overrides.append(c)

# Tune the remaining cases
for ov in overrides:
    ov["name"] = work_name(base["prefix"], ov)
    ov["job_name"] = ov["name"] + "-0" 
    span_size = 16
    nelm = ov["nodes"] * 2048
    aspect = nelm*4 / span_size**3
    ov["shape_mesh"] = [span_size, span_size, int(aspect / 4 * span_size)]
    ov["extent_mesh"] = [0.5, 0.5, 1.0 * aspect]
    ov["procs"] = 64*ov["nodes"]
    ov["io_files"] = -ov["nodes"]

from os.path import join
results = []
for ov in overrides:
    workdir = join(getcwd(), ov["name"])
    config = configure(base, ov, workdir)
    input_files = prepare(config, tusr)
    res = run(config, nekmpi_path, dep=input_files)
    results.append(res)

final = run_all(results, base)
