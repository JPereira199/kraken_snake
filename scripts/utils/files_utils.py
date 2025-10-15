from __future__ import annotations
from pathlib import Path
from typing import Iterable, Union, List
import pandas as pd
import glob

PathLike = Union[str, Path, Iterable[Union[str, Path]]]

def _coerce_to_paths(src: PathLike) -> List[Path]:
    """Accept str/Path, list/tuple of paths, or globs; return a flat list of Path objects."""
    if isinstance(src, (str, Path)):
        src = [src]
    elif not isinstance(src, Iterable):
        raise TypeError("Expected a path, glob, or an iterable of them.")

    paths: List[Path] = []
    for item in src:
        p = str(item)
        # expand globs like "*.tsv"
        matches = glob.glob(p)
        if matches:
            paths.extend(Path(m) for m in matches)
        else:
            paths.append(Path(p))
    if not paths:
        raise FileNotFoundError("No files found for the provided input(s).")
    return paths

def load_counts(paths: PathLike, sep: str = "\t", index_col: int = 0, transpose: bool = False) -> pd.DataFrame:
    """
    Load one or many count matrices (features in rows), transpose each (samples as rows → columns),
    outer-join across files on the index, and fill missing with 0.
    """
    files = _coerce_to_paths(paths)
    if transpose:
        frames = [pd.read_csv(f, sep=sep, index_col=index_col) for f in files]    
    else:
        frames = [pd.read_csv(f, sep=sep, index_col=index_col).T for f in files]
        
    df = pd.concat(frames, axis=1, join="outer")
    # Ensure numeric where possible; then fill NAs with 0
    df = df.apply(pd.to_numeric, errors="ignore").fillna(0)
    return df

def load_metadata(paths: PathLike, sep: str = "\t", id_col: str = "work_id") -> pd.DataFrame:
    """
    Load one or many metadata tables, set index to `id_col`, row-bind (axis=0),
    and drop duplicate rows by index (keep the last occurrence).
    """
    files = _coerce_to_paths(paths)
    frames = [pd.read_csv(f, sep=sep).set_index(id_col) for f in files]
    meta = pd.concat(frames, axis=0, join="outer")
    # If the same sample appears multiple times, keep the last (change to 'first' if you prefer)
    meta = meta[~meta.index.duplicated(keep="last")]
    return meta
