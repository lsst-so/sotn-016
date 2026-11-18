"""Report on outliers in the residuals between the TMA elevation and the
M1M3 inclinometer angle over a range of ``dayObs``.

Topics used, as defined in the ts-xml SAL interfaces
(https://ts-xml.lsst.io/sal_interfaces/):

- ``MTMount_elevation`` (telemetry) -- field ``actualPosition`` [deg]
- ``MTM1M3_inclinometerData`` (telemetry) -- field ``inclinometerAngle`` [deg]
"""

import argparse
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lsst.summit.utils.efdUtils import (
    getEfdData,
    makeEfdClient,
)

from lsst.summit.utils.dateTime import (
    calcNextDay,
    getDayObsEndTime,
    getDayObsStartTime,
)

__all__ = [
    "get_day_range",
    "compute_residuals",
    "flag_outliers",
    "cluster_outliers",
    "build_day_report",
    "plot_day_residuals",
    "generate_report",
]

ELEVATION_TOPIC = "lsst.sal.MTMount.elevation"
ELEVATION_FIELD = "actualPosition"
INCLINOMETER_TOPIC = "lsst.sal.MTM1M3.inclinometerData"
INCLINOMETER_FIELD = "inclinometerAngle"


def get_day_range(day_start, day_end):
    """Build the list of ``dayObs`` between ``day_start`` and ``day_end``.

    Parameters
    ----------
    day_start : int
        First ``dayObs`` of the range, in ``YYYYMMDD`` format.
    day_end : int
        Last ``dayObs`` of the range, in ``YYYYMMDD`` format (inclusive).

    Returns
    -------
    day_range : list of int
        The days in the specified date range.
    """
    if day_start > day_end:
        raise ValueError("day_start is after day_end")

    day_range = []
    day = day_start
    while day <= day_end:
        day_range.append(day)
        day = calcNextDay(day)

    return day_range


def compute_residuals(client, day_obs, tolerance="1s"):
    """Query TMA elevation and M1M3 inclinometer angle for a ``dayObs`` and
    compute their residual.

    The two telemetry streams are published at different cadences, so each
    inclinometer sample is paired with the nearest-in-time elevation sample,
    within ``tolerance``.

    Parameters
    ----------
    client : `lsst_efd_client.EfdClient`
        EFD client used to query the data.
    day_obs : int
        Day to query, in ``YYYYMMDD`` format.
    tolerance : str, optional
        Maximum time difference allowed when pairing elevation and
        inclinometer samples, as a pandas offset alias (default: ``"1s"``).

    Returns
    -------
    merged : `pandas.DataFrame`
        Paired samples indexed by time, with columns ``actualPosition``,
        ``inclinometerAngle``, and ``residual`` (``actualPosition`` minus
        ``inclinometerAngle``, in degrees). Empty if no data is available
        for either topic.
    """
    start_time = getDayObsStartTime(day_obs)
    end_time = getDayObsEndTime(day_obs)

    elevation_df = getEfdData(
        client=client,
        topic=ELEVATION_TOPIC,
        columns=[ELEVATION_FIELD],
        begin=start_time,
        end=end_time,
    )
    inclinometer_df = getEfdData(
        client=client,
        topic=INCLINOMETER_TOPIC,
        columns=[INCLINOMETER_FIELD],
        begin=start_time,
        end=end_time,
    )

    if elevation_df.empty or inclinometer_df.empty:
        return pd.DataFrame(columns=[ELEVATION_FIELD, INCLINOMETER_FIELD, "residual"])

    merged = pd.merge_asof(
        inclinometer_df.sort_index(),
        elevation_df.sort_index(),
        left_index=True,
        right_index=True,
        direction="nearest",
        tolerance=pd.Timedelta(tolerance),
    ).dropna(subset=[ELEVATION_FIELD, INCLINOMETER_FIELD])

    merged["residual"] = merged[ELEVATION_FIELD] - merged[INCLINOMETER_FIELD]

    return merged


def flag_outliers(residuals, n_sigma=5.0):
    """Flag significant outliers in a series of residuals using a robust
    (median / MAD-based) z-score.

    Parameters
    ----------
    residuals : array-like
        Residual values.
    n_sigma : float, optional
        Number of robust standard deviations away from the median beyond
        which a sample is flagged as an outlier (default: 5.0).

    Returns
    -------
    is_outlier : `numpy.ndarray` of bool
        Mask of the same length as ``residuals``, True where flagged.
    z_score : `numpy.ndarray` of float
        Robust z-score for each sample.
    """
    residuals = np.asarray(residuals, dtype=float)

    median = np.median(residuals)
    mad = np.median(np.abs(residuals - median))
    robust_sigma = 1.4826 * mad

    if robust_sigma == 0:
        robust_sigma = np.std(residuals)

    if robust_sigma == 0:
        return np.zeros(len(residuals), dtype=bool), np.zeros(len(residuals))

    z_score = (residuals - median) / robust_sigma
    is_outlier = np.abs(z_score) > n_sigma

    return is_outlier, z_score


def cluster_outliers(times, is_outlier, max_gap="60s"):
    """Group consecutive flagged samples into contiguous outlier events.

    Parameters
    ----------
    times : `pandas.DatetimeIndex`
        Timestamps of the samples, sorted ascending.
    is_outlier : array-like of bool
        Outlier mask, same length as ``times``.
    max_gap : str, optional
        Maximum gap between two flagged samples for them to be considered
        part of the same event, as a pandas offset alias (default: ``"60s"``).

    Returns
    -------
    events : list of dict
        One entry per outlier event, with keys ``start``, ``end``,
        ``n_samples``.
    """
    max_gap = pd.Timedelta(max_gap)
    is_outlier = np.asarray(is_outlier, dtype=bool)
    outlier_times = pd.DatetimeIndex(times)[is_outlier]

    if len(outlier_times) == 0:
        return []

    events = []
    event_start = outlier_times[0]
    event_end = outlier_times[0]
    n_samples = 1

    for t in outlier_times[1:]:
        if (t - event_end) <= max_gap:
            event_end = t
            n_samples += 1
        else:
            events.append({"start": event_start, "end": event_end, "n_samples": n_samples})
            event_start = t
            event_end = t
            n_samples = 1

    events.append({"start": event_start, "end": event_end, "n_samples": n_samples})

    return events


def build_day_report(
    day_obs, merged, is_outlier, z_score, max_gap="60s", min_duration=0.0, min_residual=0.0
):
    """Build a per-event outlier report for a single ``dayObs``.

    Parameters
    ----------
    day_obs : int
        Day the data belongs to, in ``YYYYMMDD`` format.
    merged : `pandas.DataFrame`
        Output of `compute_residuals` for this day.
    is_outlier : array-like of bool
        Outlier mask aligned with ``merged``.
    z_score : array-like of float
        Robust z-score aligned with ``merged``.
    max_gap : str, optional
        Passed through to `cluster_outliers` (default: ``"60s"``).
    min_duration : float, optional
        Minimum event duration, in seconds, for it to be kept in the report
        (default: 0.0, i.e. keep all events).
    min_residual : float, optional
        Minimum absolute residual, in degrees, that an event must reach
        (i.e. its ``max_abs_residual_deg``) for it to be kept in the report
        (default: 0.0, i.e. keep all events).

    Returns
    -------
    report : `pandas.DataFrame`
        One row per outlier event, with columns ``day_obs``, ``start``,
        ``end``, ``duration_s``, ``n_samples``, ``max_abs_residual_deg``,
        ``mean_residual_deg``, ``max_abs_z_score``. Empty if no events
        survive ``min_duration`` and ``min_residual``.
    """
    events = cluster_outliers(merged.index, is_outlier, max_gap=max_gap)

    rows = []
    for event in events:
        duration_s = (event["end"] - event["start"]).total_seconds()
        if duration_s < min_duration:
            continue

        window = merged.loc[event["start"] : event["end"]]
        window_z = np.asarray(z_score)[
            (merged.index >= event["start"]) & (merged.index <= event["end"])
        ]

        max_abs_residual = window["residual"].abs().max()
        if max_abs_residual < min_residual:
            continue

        rows.append(
            {
                "day_obs": day_obs,
                "start": event["start"],
                "end": event["end"],
                "duration_s": duration_s,
                "n_samples": event["n_samples"],
                "max_abs_residual_deg": max_abs_residual,
                "mean_residual_deg": window["residual"].mean(),
                "max_abs_z_score": np.abs(window_z).max(),
            }
        )

    return pd.DataFrame(rows)


def plot_day_residuals(day_obs, merged, is_outlier, outdir):
    """Plot the elevation/inclinometer residual for a ``dayObs``, with
    outlier samples highlighted, and save the figure.

    Parameters
    ----------
    day_obs : int
        Day the data belongs to, in ``YYYYMMDD`` format.
    merged : `pandas.DataFrame`
        Output of `compute_residuals` for this day.
    is_outlier : array-like of bool
        Outlier mask aligned with ``merged``.
    outdir : `pathlib.Path`
        Directory the figure is saved to.

    Returns
    -------
    plot_path : `pathlib.Path`
        Path of the saved figure.
    """
    is_outlier = np.asarray(is_outlier, dtype=bool)

    fig, ax = plt.subplots(figsize=(11, 5))

    ax.plot(merged.index, merged["residual"], color="tab:blue", lw=0.8, label="residual")
    ax.scatter(
        merged.index[is_outlier],
        merged["residual"][is_outlier],
        color="tab:red",
        s=15,
        zorder=3,
        label="outlier",
    )
    ax.axhline(0, color="black", lw=0.8, linestyle=":")

    ax.set_xlabel("Time (UTC)")
    ax.set_ylabel("TMA elevation - M1M3 inclinometer angle (deg)")
    ax.set_title(f"TMA elevation / M1M3 inclinometer residual - dayObs {day_obs}")
    ax.legend()
    fig.autofmt_xdate()
    fig.tight_layout()

    plot_path = Path(outdir) / f"residuals_{day_obs}.png"
    fig.savefig(plot_path)
    plt.close(fig)

    return plot_path


def generate_report(
    day_start,
    day_end,
    n_sigma=5.0,
    tolerance="1s",
    max_gap="60s",
    min_duration=0.0,
    min_residual=0.0,
    outdir=None,
    client=None,
):
    """Query TMA elevation and M1M3 inclinometer data for a range of
    ``dayObs``, flag significant outliers in their residual, and produce a
    report and plots for the days where outliers were found.

    Parameters
    ----------
    day_start : int
        First ``dayObs`` of the range, in ``YYYYMMDD`` format.
    day_end : int
        Last ``dayObs`` of the range, in ``YYYYMMDD`` format (inclusive).
    n_sigma : float, optional
        Outlier threshold passed to `flag_outliers` (default: 5.0).
    tolerance : str, optional
        Passed through to `compute_residuals` (default: ``"1s"``).
    max_gap : str, optional
        Passed through to `cluster_outliers` (default: ``"60s"``).
    min_duration : float, optional
        Passed through to `build_day_report` (default: 0.0).
    min_residual : float, optional
        Passed through to `build_day_report` (default: 0.0).
    outdir : str or `pathlib.Path`, optional
        Directory the report and plots are written to. If not provided, a
        timestamped directory is created in the current working directory.
    client : `lsst_efd_client.EfdClient`, optional
        EFD client to use. If not provided, one is created with
        `makeEfdClient`.

    Returns
    -------
    report : `pandas.DataFrame`
        One row per outlier event found in the range, sorted by start time.
        Empty if no outliers were found.
    outdir : `pathlib.Path`
        Directory the report CSV and the per-day residual plots were
        written to.
    """
    if client is None:
        client = makeEfdClient()

    if outdir is None:
        timestamp = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        outdir = f"tma_elevation_inclinometer_report_{timestamp}"
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    day_range = get_day_range(day_start, day_end)

    day_reports = []
    for day_obs in day_range:
        merged = compute_residuals(client, day_obs, tolerance=tolerance)

        if merged.empty:
            print(f"dayObs {day_obs}: no overlapping elevation/inclinometer data, skipping")
            continue

        is_outlier, z_score = flag_outliers(merged["residual"], n_sigma=n_sigma)
        day_report = build_day_report(
            day_obs,
            merged,
            is_outlier,
            z_score,
            max_gap=max_gap,
            min_duration=min_duration,
            min_residual=min_residual,
        )

        if day_report.empty:
            print(f"dayObs {day_obs}: {len(merged)} samples, no significant outliers")
            continue

        plot_path = plot_day_residuals(day_obs, merged, is_outlier, outdir)
        print(
            f"dayObs {day_obs}: {len(merged)} samples, "
            f"{len(day_report)} outlier event(s), plot saved to {plot_path}"
        )
        day_reports.append(day_report)

    if day_reports:
        report = pd.concat(day_reports, ignore_index=True).sort_values("start").reset_index(drop=True)
    else:
        report = pd.DataFrame(
            columns=[
                "day_obs",
                "start",
                "end",
                "duration_s",
                "n_samples",
                "max_abs_residual_deg",
                "mean_residual_deg",
                "max_abs_z_score",
            ]
        )

    report_path = outdir / "outlier_report.csv"
    report.to_csv(report_path, index=False)

    print(f"\n{len(report)} outlier event(s) found between dayObs {day_start} and {day_end}")
    print(f"Report written to {report_path}")
    if not report.empty:
        print(report.to_string(index=False))

    return report, outdir


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Report significant outliers in the residual between TMA elevation "
            "and M1M3 inclinometer angle over a range of dayObs, and plot the "
            "residuals for the days where outliers were found."
        )
    )
    parser.add_argument(
        "--start-day",
        type=int,
        required=True,
        help="First dayObs of the range, in YYYYMMDD format",
    )
    parser.add_argument(
        "--end-day",
        type=int,
        required=True,
        help="Last dayObs of the range (inclusive), in YYYYMMDD format",
    )
    parser.add_argument(
        "--n-sigma",
        type=float,
        default=5.0,
        help="Robust (median/MAD) sigma threshold for flagging outliers (default: 5.0)",
    )
    parser.add_argument(
        "--tolerance",
        type=str,
        default="1s",
        help="Max time difference when pairing elevation/inclinometer samples (default: 1s)",
    )
    parser.add_argument(
        "--max-gap",
        type=str,
        default="60s",
        help="Max gap between flagged samples to group them into one event (default: 60s)",
    )
    parser.add_argument(
        "--min-duration",
        type=float,
        default=0.0,
        help="Minimum outlier event duration, in seconds, to include in the report (default: 0.0)",
    )
    parser.add_argument(
        "--min-residual",
        type=float,
        default=0.0,
        help=(
            "Minimum absolute residual, in degrees, an outlier event must reach "
            "to be included in the report (default: 0.0)"
        ),
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=None,
        help="Output directory for the report and plots (default: timestamped directory)",
    )
    args = parser.parse_args()

    generate_report(
        args.start_day,
        args.end_day,
        n_sigma=args.n_sigma,
        tolerance=args.tolerance,
        max_gap=args.max_gap,
        min_duration=args.min_duration,
        min_residual=args.min_residual,
        outdir=args.outdir,
    )


if __name__ == "__main__":
    main()
