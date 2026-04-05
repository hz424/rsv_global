import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


TREE_CONFIG = {
    "A": {"json": "treeA/rsv_A_3.json"},
    "B": {"json": "treeB/rsv_B_5.json"},
}


# These are the only country/state labels currently observed in the local trees.

COUNTRY_NORMALIZATION: Dict[str, Dict[str, str]] = {
    "Australia": {"display": "Australia", "map": "Australia"},
    "Beijing": {"display": "Beijing", "map": "China"},
    "Canada": {"display": "Canada", "map": "Canada"},
    "China": {"display": "China", "map": "China"},
    "Cote_dIvoire": {"display": "Cote d'Ivoire", "map": "Ivory Coast"},
    "England": {"display": "England", "map": "United Kingdom"},
    "France": {"display": "France", "map": "France"},
    "India": {"display": "India", "map": "India"},
    "Iran": {"display": "Iran", "map": "Iran"},
    "Ireland": {"display": "Ireland", "map": "Ireland"},
    "Italy": {"display": "Italy", "map": "Italy"},
    "Kenya": {"display": "Kenya", "map": "Kenya"},
    "Norway": {"display": "Norway", "map": "Norway"},
    "Pakistan": {"display": "Pakistan", "map": "Pakistan"},
    "Philippines": {"display": "Philippines", "map": "Philippines"},
    "Qatar": {"display": "Qatar", "map": "Qatar"},
    "Russia": {"display": "Russia", "map": "Russia"},
    "South_Africa": {"display": "South Africa", "map": "South Africa"},
    "Spain": {"display": "Spain", "map": "Spain"},
    "Thailand": {"display": "Thailand", "map": "Thailand"},
    "USA": {"display": "USA", "map": "United States"},
}


PANEL_ORDER = [
    ("A", "Import", "Introductions RSV A"),
    ("B", "Import", "Introductions RSV B"),
    ("A", "Export", "Exports RSV A"),
    ("B", "Export", "Exports RSV B"),
]


def display_country(name: str) -> str:
    normalized = COUNTRY_NORMALIZATION.get(name)
    if normalized:
        return normalized["display"]
    return name.replace("_", " ")


def map_country_name(name: str) -> str:
    normalized = COUNTRY_NORMALIZATION.get(name)
    if normalized:
        return normalized["map"]
    return display_country(name)


def validate_country_normalization(edges: pd.DataFrame) -> None:
    observed = set(edges["partner_country"].dropna().astype(str))
    missing = sorted(country for country in observed if country not in COUNTRY_NORMALIZATION)
    if missing:
        raise ValueError(
            "Unmapped partner_country values detected. "
            "Extend COUNTRY_NORMALIZATION before plotting: %s" % ", ".join(missing)
        )


def get_country_state(node: dict) -> Tuple[Optional[str], Optional[float]]:
    country = node.get("node_attrs", {}).get("country", {})
    if not isinstance(country, dict):
        return None, None
    value = country.get("value")
    confidence = country.get("confidence", {})
    posterior = confidence.get(value) if isinstance(confidence, dict) else None
    return value, posterior


def get_num_date(node: dict) -> Optional[float]:
    nd = node.get("node_attrs", {}).get("num_date", {})
    if isinstance(nd, dict):
        return nd.get("value")
    return None


def extract_uae_transition_edges(subtype: str, json_path: Path, threshold: float) -> pd.DataFrame:
    tree = json.loads(json_path.read_text())["tree"]
    rows = []  # type: List[dict]
    stack = [tree]

    while stack:
        parent = stack.pop()
        parent_country, parent_posterior = get_country_state(parent)
        parent_date = get_num_date(parent)

        for child in parent.get("children", []):
            child_country, child_posterior = get_country_state(child)
            child_date = get_num_date(child)

            if not parent_country or not child_country or parent_country == child_country:
                stack.append(child)
                continue

            if "UAE" not in (parent_country, child_country):
                stack.append(child)
                continue

            direction = "Import" if child_country == "UAE" else "Export"
            partner_country = parent_country if direction == "Import" else child_country
            keep = (
                parent_posterior is not None
                and child_posterior is not None
                and parent_posterior >= threshold
                and child_posterior >= threshold
            )

            rows.append(
                {
                    "subtype": subtype,
                    "direction": direction,
                    "partner_country": partner_country,
                    "partner_country_label": display_country(partner_country),
                    "map_country": map_country_name(partner_country),
                    "parent_country": parent_country,
                    "child_country": child_country,
                    "parent_posterior": parent_posterior,
                    "child_posterior": child_posterior,
                    "parent_num_date": parent_date,
                    "child_num_date": child_date,
                    "parent_node": parent.get("name"),
                    "child_node": child.get("name"),
                    "posterior_threshold": threshold,
                    "kept_strict": keep,
                }
            )

            stack.append(child)

    return pd.DataFrame(rows)


def build_country_summary(edges: pd.DataFrame) -> pd.DataFrame:
    loose = (
        edges.groupby(
            ["subtype", "direction", "partner_country", "partner_country_label", "map_country"],
            as_index=False,
        )
        .size()
        .rename(columns={"size": "count_loose"})
    )
    strict = (
        edges.loc[edges["kept_strict"]]
        .groupby(
            ["subtype", "direction", "partner_country", "partner_country_label", "map_country"],
            as_index=False,
        )
        .size()
        .rename(columns={"size": "count_strict"})
    )
    summary = loose.merge(
        strict,
        on=["subtype", "direction", "partner_country", "partner_country_label", "map_country"],
        how="outer",
    ).fillna({"count_loose": 0, "count_strict": 0})
    summary["count_loose"] = summary["count_loose"].astype(int)
    summary["count_strict"] = summary["count_strict"].astype(int)
    summary["filtered_out"] = summary["count_loose"] - summary["count_strict"]
    return summary.sort_values(["subtype", "direction", "count_strict", "partner_country"])


def build_map_summary(summary: pd.DataFrame) -> pd.DataFrame:
    strict = summary.loc[summary["count_strict"] > 0].copy()
    if strict.empty:
        return strict

    grouped = (
        strict.groupby(["subtype", "direction", "map_country"], as_index=False)["count_strict"]
        .sum()
        .rename(columns={"count_strict": "map_count"})
    )
    grouped["map_count_label"] = grouped["map_count"].astype(int).astype(str)
    return grouped.sort_values(["subtype", "direction", "map_count", "map_country"])


def add_map_trace(fig, panel_df: pd.DataFrame, row: int, col: int, colorbar_x: float) -> None:
    zmax = max(panel_df["map_count"].max(), 1) if not panel_df.empty else 1
    trace = go.Choropleth(
        locations=panel_df["map_country"] if not panel_df.empty else [],
        locationmode="country names",
        z=panel_df["map_count"] if not panel_df.empty else [],
        text=panel_df["map_country"] if not panel_df.empty else [],
        customdata=panel_df[["map_count_label"]] if not panel_df.empty else None,
        colorscale="Viridis",
        zmin=1,
        zmax=zmax,
        marker_line_color="#4f4f4f",
        marker_line_width=0.6,
        colorbar=dict(
            title="count",
            len=0.35,
            thickness=12,
            x=colorbar_x,
            y=0.78 if row == 1 else 0.25,
            xanchor="left",
            outlinewidth=0,
        ),
        hovertemplate="%{text}<br>count=%{z}<extra></extra>",
        showscale=True,
    )
    fig.add_trace(trace, row=row, col=col)


def format_geo(fig) -> None:
    fig.update_geos(
        projection_type="equirectangular",
        projection_scale=1.0,
        showcountries=True,
        countrycolor="#5b5b5b",
        showcoastlines=False,
        showland=True,
        landcolor="#e8f0fb",
        showocean=True,
        oceancolor="white",
        showframe=True,
        framecolor="#5b5b5b",
        framewidth=1.0,
        lataxis_showgrid=False,
        lataxis_range=[-60, 88],
        lonaxis_showgrid=False,
        lonaxis_range=[-180, 180],
    )


def make_global_map(map_summary: pd.DataFrame, threshold: float, out_base: Path) -> None:
    fig = make_subplots(
        rows=2,
        cols=2,
        specs=[[{"type": "choropleth"}, {"type": "choropleth"}], [{"type": "choropleth"}, {"type": "choropleth"}]],
        subplot_titles=[title for _, _, title in PANEL_ORDER],
        horizontal_spacing=0.015,
        vertical_spacing=0.08,
    )

    colorbar_x = {(1, 1): 0.455, (1, 2): 0.952, (2, 1): 0.455, (2, 2): 0.952}

    for idx, (subtype, direction, _) in enumerate(PANEL_ORDER, start=1):
        row = 1 if idx <= 2 else 2
        col = 1 if idx in (1, 3) else 2
        panel = map_summary.loc[(map_summary["subtype"] == subtype) & (map_summary["direction"] == direction)]
        add_map_trace(fig, panel, row=row, col=col, colorbar_x=colorbar_x[(row, col)])

    format_geo(fig)

    fig.update_layout(
        geo=dict(domain={"x": [0.00, 0.448], "y": [0.56, 0.98]}),
        geo2=dict(domain={"x": [0.505, 0.95], "y": [0.56, 0.98]}),
        geo3=dict(domain={"x": [0.00, 0.448], "y": [0.07, 0.49]}),
        geo4=dict(domain={"x": [0.505, 0.95], "y": [0.07, 0.49]}),
    )

    fig.update_layout(
        width=1500,
        height=690,
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin=dict(l=8, r=8, t=36, b=44),
        font=dict(family="Arial, Liberation Sans, DejaVu Sans", size=16, color="#3f5679"),
        annotations=[
            dict(
                x=0.5,
                y=-0.01,
                xref="paper",
                yref="paper",
                showarrow=False,
                text=(
                    "Discrete parent→child country-state transition edges involving the UAE. "
                    "Parent and child assigned-country posterior both ≥ %.1f."
                )
                % threshold,
                font=dict(size=15, color="#3f5679"),
            )
        ]
        + list(fig.layout.annotations),
    )

    fig.write_image(str(out_base.with_suffix(".png")), scale=2)
    fig.write_image(str(out_base.with_suffix(".svg")))
    fig.write_image(str(out_base.with_suffix(".pdf")))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Make posterior-filtered UAE import/export world maps.")
    parser.add_argument("--threshold", type=float, default=0.9)
    parser.add_argument("--tree-a-json", type=Path, default=None)
    parser.add_argument("--tree-b-json", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[2]
    workdir = Path(__file__).resolve().parents[1]
    figure_dir = workdir / "figures"
    table_dir = workdir / "tables"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    edge_tables = []
    for subtype, cfg in TREE_CONFIG.items():
        if subtype == "A" and args.tree_a_json is not None:
            json_path = args.tree_a_json
        elif subtype == "B" and args.tree_b_json is not None:
            json_path = args.tree_b_json
        else:
            json_path = root / cfg["json"]
        edge_tables.append(
            extract_uae_transition_edges(
                subtype=subtype,
                json_path=json_path,
                threshold=args.threshold,
            )
        )

    edges = pd.concat(edge_tables, ignore_index=True)
    validate_country_normalization(edges)
    summary = build_country_summary(edges)
    map_summary = build_map_summary(summary)

    edges.to_csv(table_dir / "uae_transition_edges.tsv", sep="\t", index=False)
    summary.to_csv(table_dir / "uae_transition_country_counts.tsv", sep="\t", index=False)
    map_summary.to_csv(table_dir / "uae_transition_country_counts_mapped.tsv", sep="\t", index=False)

    make_global_map(
        map_summary=map_summary,
        threshold=args.threshold,
        out_base=figure_dir / "Fig5_posterior_filtered_global_map",
    )


if __name__ == "__main__":
    main()
