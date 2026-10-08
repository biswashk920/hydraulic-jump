"""Reading channel descriptions (JSON or CSV) and writing results."""
import json
import os
import pandas as pd
from .sections import ChannelError


def parse_bc(text):
    """'depth:2.5', 'gate:0.4,0.6', 'weir:1.2,1.7', 'level:101.5', 'normal', 'critical', 'overfall', 'none'."""
    t, _, vals = text.partition(":")
    t = t.strip().lower()
    v = [s for s in vals.split(",") if s.strip()]
    try:
        v = [float(s) for s in v]
    except ValueError:
        raise ChannelError(f"could not read the numbers in '{text}'.")
    names = {"depth": ["depth"], "level": ["level"], "gate": ["opening", "cc"], "weir": ["crest", "cw"]}
    bc = {"type": t}
    for k, val in zip(names.get(t, []), v):
        bc[k] = val
    return bc


def _csv_to_spec(path):
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    need = {"length", "slope", "n", "section"}
    if not need <= set(df.columns):
        raise ChannelError(f"the CSV needs the columns: length, slope, n, section (and width, side_slope, diameter "
                           f"as needed). Found: {', '.join(df.columns)}.")
    reaches = []
    for _, row in df.iterrows():
        sec = {"type": str(row["section"]).strip().lower()}
        for k in ("width", "side_slope", "diameter"):
            if k in df.columns and pd.notna(row[k]):
                sec[k] = float(row[k])
        reaches.append({"length": float(row["length"]), "slope": float(row["slope"]), "n": float(row["n"]), "section": sec})
    return {"name": os.path.splitext(os.path.basename(path))[0], "reaches": reaches}


def load_spec(path):
    if not os.path.exists(path):
        raise ChannelError(f"file not found: {path}")
    if path.lower().endswith(".csv"):
        return _csv_to_spec(path)
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        raise ChannelError(f"{path} is not valid JSON (line {e.lineno}, column {e.colno}): {e.msg}.")
