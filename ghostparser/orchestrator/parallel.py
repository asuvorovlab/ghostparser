"""Worker-pool plumbing shared by every parallel stage: the CPU count, the
start-method choice and an order-preserving map over shared read-only state.
"""

import multiprocessing as mp
import os
import time
from multiprocessing import cpu_count

# Worker-global state, seeded by the pool initializer or inherited via fork.
_WORKER_FUNC = None
_WORKER_SHARED = None


def available_cpu_count():
    """Count the CPUs this process may run on, from its affinity mask rather
    than the machine's core count.

    Returns:
        The usable CPU count, at least 1.
    """
    counter = getattr(os, "process_cpu_count", None)
    if counter is not None:
        return max(1, counter() or 1)
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, cpu_count())


def resolve_worker_count(total_items, processes):
    """Resolve the worker count from the requested processes and item count.

    Args:
        total_items: Number of work items available.
        processes: Requested worker count; ``0`` means every CPU available to
            this process.

    Returns:
        The clamped worker count (at least 1, at most ``total_items``).
    """
    worker_count = processes or available_cpu_count()
    return max(1, min(worker_count, total_items))


def get_mp_context(prefer_fork=None):
    """Select a multiprocessing context, preferring fork then forkserver then spawn.

    Args:
        prefer_fork: Whether to prefer ``fork`` when available. Defaults to
            ``True`` on POSIX platforms when ``None``.

    Returns:
        A multiprocessing context (or the ``multiprocessing`` module itself).
    """
    if prefer_fork is None:
        prefer_fork = os.name == "posix"

    if hasattr(mp, "get_context"):
        methods = mp.get_all_start_methods()
        if prefer_fork and "fork" in methods:
            return mp.get_context("fork")
        if "forkserver" in methods:
            return mp.get_context("forkserver")
        if "spawn" in methods:
            return mp.get_context("spawn")
    return mp


def _init_worker(func, shared, inherited):
    """Seed a worker's globals with the mapped function and the shared state.

    Args:
        func: The function every item is passed to.
        shared: The read-only state handed to ``func``; ignored when
            ``inherited``.
        inherited: Whether the worker already holds the state via fork.
    """
    global _WORKER_FUNC, _WORKER_SHARED
    _WORKER_FUNC = func
    if not inherited:
        _WORKER_SHARED = shared


def _run_item(item):
    """Apply the worker's function to one item, timing its CPU.

    Args:
        item: One work item.

    Returns:
        A tuple ``(result, cpu_seconds)``.
    """
    start_cpu = time.process_time()
    result = _WORKER_FUNC(_WORKER_SHARED, item)
    return result, time.process_time() - start_cpu


def map_ordered(func, items, worker_count, shared=None, chunksize=None):
    """Apply ``func(shared, item)`` to every item, in input order, across workers.

    Under fork the workers inherit ``shared`` copy-on-write instead of
    receiving a pickled copy each.

    Args:
        func: A module-level function taking ``(shared, item)``.
        items: The work items.
        worker_count: Number of worker processes; ``1`` or fewer runs serially
            in this process.
        shared: Read-only state every call receives.
        chunksize: Items per dispatched task; defaults to about four tasks per
            worker.

    Returns:
        A tuple ``(results, worker_cpu_seconds)``: the results in input order,
        and the CPU time spent in workers, which the caller's own
        ``process_time`` cannot see (``0.0`` when run serially).
    """
    global _WORKER_SHARED
    items = list(items)
    if worker_count <= 1 or len(items) <= 1:
        return [func(shared, item) for item in items], 0.0

    if chunksize is None:
        chunksize = max(1, len(items) // (worker_count * 4))
    ctx = get_mp_context()
    inherited = (
        hasattr(ctx, "get_start_method") and ctx.get_start_method() == "fork"
    )
    if inherited:
        _WORKER_SHARED = shared
    try:
        with ctx.Pool(
            processes=worker_count,
            initializer=_init_worker,
            initargs=(func, None if inherited else shared, inherited),
        ) as pool:
            payloads = list(pool.imap(_run_item, items, chunksize=chunksize))
    finally:
        if inherited:
            _WORKER_SHARED = None
    return [result for result, _ in payloads], sum(cpu for _, cpu in payloads)
