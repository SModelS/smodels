#!/usr/bin/env python3

from __future__ import annotations

import argparse
import itertools
import random
import time
from collections import defaultdict

from smodels.base import runtime
from smodels.base.model import Model
from smodels.base.physicsUnits import GeV, fb
from smodels.base.smodelsLogging import setLogLevel
from smodels.decomposition.decomposer import decompose
from smodels.decomposition.decomposerNew import decomposeNew
from smodels.decomposition.decomposerOld import decomposeOld
from smodels.share.models.SMparticles import SMList
from smodels.tools.particlesLoader import load


RUNS = [
    ("Old (BFS)", decomposeOld),
    ("Current (DFS)", decompose),
    ("New (subtree cache)", decomposeNew),
]


def clear_particle_caches(model: Model) -> None:
    for particle in model.BSMparticles:
        for attr in ("_decayTrees", "_maxDecayBR"):
            try:
                delattr(particle, attr)
            except AttributeError:
                pass


def load_model(runtime_model_file: str, input_file: str) -> Model:
    runtime.modelFile = runtime_model_file
    bsm_list = load()
    model = Model(BSMparticles=bsm_list, SMparticles=SMList)
    model.updateParticles(
        inputFile=input_file,
        ignorePromptQNumbers=["eCharge", "spin"],
    )
    return model


def benchmark_decomposers(
    model_factory,
    sigmacut,
    cycles: int,
    clear_before_each_call: bool,
    seed: int,
    *,
    fresh_model_per_call: bool,
    **kwargs,
):
    rng = random.Random(seed)
    results = defaultdict(list)
    orders = list(itertools.permutations(RUNS)) * cycles
    rng.shuffle(orders)
    shared_model = None if fresh_model_per_call else model_factory()

    for irep, order in enumerate(orders, start=1):
        labels = [label for label, _ in order]
        print(f"rep {irep}/{len(orders)} order = {labels}")
        for position, (label, fn) in enumerate(order, start=1):
            model = model_factory() if fresh_model_per_call else shared_model
            if clear_before_each_call:
                clear_particle_caches(model)
            t0 = time.perf_counter()
            result = fn(model, sigmacut, **kwargs)
            dt = time.perf_counter() - t0
            n_topologies = len(result)
            n_sms = len(result.getSMSList())
            results[label].append(
                {
                    "time": dt,
                    "position": position,
                    "topologies": n_topologies,
                    "sms": n_sms,
                }
            )
            print(
                f"  {label:<20} time={dt:8.3f}s topologies={n_topologies:3d} "
                f"sms={n_sms:5d} pos={position}"
            )

    return results


def _mean(values) -> float:
    return sum(values) / len(values)


def _median(values) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return 0.5 * (ordered[middle - 1] + ordered[middle])


def summarize_benchmark(results, title: str) -> None:
    print(f"\n-- {title} --")
    print(
        f"{'Implementation':<22} {'median':>8} {'mean':>8} {'min':>8} {'max':>8} "
        f"{'pos1':>8} {'pos2':>8} {'pos3':>8} {'tops':>10} {'sms':>8}"
    )
    print("-" * 96)
    for label, _ in RUNS:
        rows = results[label]
        times = [row["time"] for row in rows]
        pos1 = [row["time"] for row in rows if row["position"] == 1]
        pos2 = [row["time"] for row in rows if row["position"] == 2]
        pos3 = [row["time"] for row in rows if row["position"] == 3]
        topologies = sorted({row["topologies"] for row in rows})
        sms_counts = sorted({row["sms"] for row in rows})
        print(
            f"{label:<22} {_median(times):8.3f} {_mean(times):8.3f} "
            f"{min(times):8.3f} {max(times):8.3f} {_mean(pos1):8.3f} "
            f"{_mean(pos2):8.3f} {_mean(pos3):8.3f} "
            f"{str(topologies):>10} {str(sms_counts):>8}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Benchmark decompose, decomposeOld, and decomposeNew across balanced call orders.",
    )
    parser.add_argument("--runtime-model-file", default="nmssmPoints/000.slha")
    parser.add_argument("--input-file", default="nmssmPoints/022.slha")
    parser.add_argument("--sigmacut-fb", type=float, default=0.05)
    parser.add_argument("--mingap-gev", type=float, default=10.0)
    parser.add_argument("--mingap-isr-gev", type=float, default=10.0)
    parser.add_argument("--cycles", type=int, default=1)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--log-level", default="info")
    parser.add_argument("--skip-warm", action="store_true")
    parser.add_argument("--skip-cold", action="store_true")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.cycles < 1:
        parser.error("--cycles must be at least 1")
    if args.skip_warm and args.skip_cold:
        parser.error("at least one of --skip-warm or --skip-cold must be false")

    setLogLevel(args.log_level)

    model_factory = lambda: load_model(args.runtime_model_file, args.input_file)
    sigmacut = args.sigmacut_fb * fb
    kwargs = {
        "massCompress": True,
        "invisibleCompress": True,
        "minmassgap": args.mingap_gev * GeV,
        "minmassgapISR": args.mingap_isr_gev * GeV,
    }

    if not args.skip_warm:
        warm_results = benchmark_decomposers(
            model_factory,
            sigmacut,
            cycles=args.cycles,
            clear_before_each_call=False,
            seed=args.seed,
            fresh_model_per_call=False,
            **kwargs,
        )
        summarize_benchmark(warm_results, "Warm-cache runs (shared model within each phase)")

    if not args.skip_cold:
        cold_results = benchmark_decomposers(
            model_factory,
            sigmacut,
            cycles=args.cycles,
            clear_before_each_call=False,
            seed=args.seed,
            fresh_model_per_call=True,
            **kwargs,
        )
        summarize_benchmark(cold_results, "Cold-start runs (fresh model for every call)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())