import glob
import html
import os
import csv
import sys
import shutil
import zipfile
import re
import datetime
import statistics

proj = sys.argv[1]
report_input_files = [
    p for p in glob.glob("**/*", recursive=True)
    if os.path.isfile(p)
]
species_model = sys.argv[2]
star_alignment_mode = sys.argv[3]
sample_barcode_sheet = sys.argv[4] if len(sys.argv) > 4 else ""
out_path = "QC_Report.html"
assets_dir = "QC_Report_assets"
per_sample_dir = "QC_Report"
os.makedirs(assets_dir, exist_ok=True)
os.makedirs(per_sample_dir, exist_ok=True)

if star_alignment_mode == "single":
    active_starsolo_roots = ["STARsolo"]
elif star_alignment_mode == "paired":
    active_starsolo_roots = ["STARsolo_paired"]
else:
    active_starsolo_roots = ["STARsolo", "STARsolo_paired"]

def hydrate_report_inputs():
    roots = ("STARsolo", "STARsolo_paired", "ATAC", "multiqc_atac")
    sample_roots = {}
    for src in report_input_files:
        norm = os.path.normpath(src)
        parts = norm.split(os.sep)
        for i, part in enumerate(parts):
            if part in ("STARsolo", "STARsolo_paired") and len(parts) > i + 1:
                sample = parts[i + 1]
                sample_roots.setdefault(sample, set()).add(part)
                break
    for src in report_input_files:
        norm = os.path.normpath(src)
        parts = norm.split(os.sep)
        rel = None
        for i, part in enumerate(parts):
            if part in roots and len(parts) > i + 1:
                rel = os.path.join(*parts[i:])
                break
        if rel is None:
            base = os.path.basename(src)
            if base.endswith("_knee_plot.png"):
                sample = base[:-len("_knee_plot.png")]
                candidate_roots = sorted(sample_roots.get(sample, set()))
                if not candidate_roots:
                    candidate_roots = list(active_starsolo_roots)
                for root in candidate_roots:
                    dst = os.path.join(proj, root, sample, base)
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    shutil.copy2(src, dst)
            continue
        dst = os.path.join(proj, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)

hydrate_report_inputs()

def ensure_knee_plots():
    work_candidates = sorted(set(
        glob.glob("**/*_knee_plot.png", recursive=True) +
        glob.glob("*_knee_plot.png")
    ))
    work_candidates = [p for p in work_candidates if os.path.isfile(p)]
    by_sample = {}
    for p in work_candidates:
        base = os.path.basename(p)
        if not base.endswith("_knee_plot.png"):
            continue
        sample = base[:-len("_knee_plot.png")]
        if sample not in by_sample:
            by_sample[sample] = p
        else:
            existing = by_sample[sample]
            if ("STARsolo" not in existing.split(os.sep)) and ("STARsolo" in p.split(os.sep)):
                by_sample[sample] = p
    for sample, src in by_sample.items():
        for root in active_starsolo_roots:
            dst_dir = os.path.join(proj, root, sample)
            os.makedirs(dst_dir, exist_ok=True)
            dst = os.path.join(dst_dir, f"{sample}_knee_plot.png")
            try:
                shutil.copy2(src, dst)
            except Exception:
                pass

ensure_knee_plots()

def ensure_atac_files():
    atac_patterns = [
        "*.q30.rmdup.flagstat.txt",
        "*.q30.rmdup.idxstats.txt",
        "*.q30.rmdup.stats.txt",
        "*.q30.mapped.flagstat.txt",
        "*.q30.mapped.idxstats.txt",
        "*.q30.mapped.stats.txt",
        "*.atac_cells.summary.tsv",
        "*.atac_cells.counts.tsv",
    ]
    name_patterns = [
        re.compile(r"(.+?)\.q30\.(?:rmdup|mapped|possort)\.(?:flagstat|idxstats|stats)\.txt\$"),
        re.compile(r"(.+?)\.atac_cells\.(?:summary|counts)\.tsv\$"),
    ]
    seen = set()
    for pat in atac_patterns:
        for src in glob.glob(f"**/{pat}", recursive=True) + glob.glob(pat):
            if not os.path.isfile(src) or src in seen:
                continue
            seen.add(src)
            base = os.path.basename(src)
            sample = None
            for rx in name_patterns:
                m = rx.match(base)
                if m:
                    sample = m.group(1)
                    break
            if not sample:
                continue
            dst_dir = os.path.join(proj, "ATAC", sample)
            os.makedirs(dst_dir, exist_ok=True)
            dst = os.path.join(dst_dir, base)
            try:
                shutil.copy2(src, dst)
            except Exception:
                pass

ensure_atac_files()

def rel_list(pattern):
    return sorted([
        os.path.relpath(p, proj)
        for p in glob.glob(os.path.join(proj, pattern))
        if os.path.isfile(p)
    ])

def rel_list_recursive(pattern):
    return sorted([
        os.path.relpath(p, proj)
        for p in glob.glob(os.path.join(proj, pattern), recursive=True)
        if os.path.isfile(p)
    ])

def safe_asset_name(rel_path):
    name = rel_path.replace("/", "__")
    return name.replace(" ", "_")

def stage_asset(rel_path):
    src = os.path.join(proj, rel_path)
    if not os.path.isfile(src):
        return None
    dest = os.path.join(assets_dir, safe_asset_name(rel_path))
    shutil.copy2(src, dest)
    return dest

def stage_asset_in(rel_path, dest_dir):
    src = os.path.join(proj, rel_path)
    if not os.path.isfile(src):
        return None
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, safe_asset_name(rel_path))
    shutil.copy2(src, dest)
    return dest

def read_table_preview(rel_path, max_rows=8):
    abs_path = os.path.join(proj, rel_path)
    if not os.path.exists(abs_path):
        return "<p><em>Missing file</em></p>"
    rows = []
    delim = "," if rel_path.lower().endswith(".csv") else chr(9)
    try:
        with open(abs_path, newline="") as fh:
            reader = csv.reader(fh, delimiter=delim)
            for i, row in enumerate(reader):
                rows.append([html.escape(x) for x in row])
                if max_rows is not None and i + 1 >= max_rows:
                    break
    except Exception as e:
        return f"<p><em>Could not parse table: {html.escape(str(e))}</em></p>"
    if not rows:
        return "<p><em>File is empty</em></p>"
    cells = []
    for ridx, row in enumerate(rows):
        tag = "th" if ridx == 0 else "td"
        cells.append("<tr>" + "".join(f"<{tag}>{c}</{tag}>" for c in row) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def load_demux_sample_names(rel_paths):
    names = set()
    delim_tsv = chr(9)
    for rel in rel_paths:
        p = os.path.join(proj, rel)
        if not os.path.isfile(p):
            continue
        delim = "," if rel.lower().endswith(".csv") else delim_tsv
        try:
            with open(p, newline="") as fh:
                reader = csv.reader(fh, delimiter=delim)
                rows = list(reader)
        except Exception:
            continue
        if not rows:
            continue
        header = rows[0]
        try:
            name_idx = header.index("Sample_Name")
        except ValueError:
            continue
        for r in rows[1:]:
            if len(r) > name_idx and r[name_idx].strip():
                names.add(r[name_idx].strip())
    return names

def demux_stats_html_for_sample(rel_paths, sample):
    delim_tsv = chr(9)
    for rel in rel_paths:
        p = os.path.join(proj, rel)
        if not os.path.isfile(p):
            continue
        delim = "," if rel.lower().endswith(".csv") else delim_tsv
        try:
            with open(p, newline="") as fh:
                reader = csv.reader(fh, delimiter=delim)
                all_rows = list(reader)
        except Exception:
            continue
        if not all_rows:
            continue
        header = all_rows[0]
        try:
            name_idx = header.index("Sample_Name")
        except ValueError:
            continue
        body = [r for r in all_rows[1:] if len(r) > name_idx and r[name_idx] == sample]
        if not body:
            continue
        esc_rows = [[html.escape(x) for x in row] for row in [header] + body]
        cells = []
        for ridx, row in enumerate(esc_rows):
            tag = "th" if ridx == 0 else "td"
            cells.append("<tr>" + "".join(f"<{tag}>{c}</{tag}>" for c in row) + "</tr>")
        return "<table>" + "".join(cells) + "</table>"
    return f"<p><em>No demultiplex stats row for <code>{html.escape(sample)}</code>.</em></p>"

def read_text_preview(rel_path, max_lines=80):
    abs_path = os.path.join(proj, rel_path)
    if not os.path.exists(abs_path):
        return "<p><em>Missing file</em></p>"
    try:
        lines = []
        with open(abs_path, "r", errors="replace") as fh:
            for i, line in enumerate(fh):
                lines.append(html.escape(line.rstrip("\n")))
                if i + 1 >= max_lines:
                    break
        if not lines:
            return "<p><em>File is empty</em></p>"
        return "<pre>" + "\n".join(lines) + "</pre>"
    except Exception as e:
        return f"<p><em>Could not read file: {html.escape(str(e))}</em></p>"

def links_block(title, paths):
    if not paths:
        return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
    items = []
    for p in paths:
        asset = stage_asset(p)
        if asset is None:
            continue
        items.append(f'<li><a href="{html.escape(asset)}">{html.escape(p)}</a></li>')
    if not items:
        return f"<h3>{html.escape(title)}</h3><p><em>No readable files found.</em></p>"
    return f"<h3>{html.escape(title)}</h3><ul>{''.join(items)}</ul>"

def text_files_block(title, paths, max_lines=80):
    if not paths:
        return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
    chunks = [f"<h3>{html.escape(title)}</h3>"]
    for p in paths:
        chunks.append(f"<h4>{html.escape(p)}</h4>")
        chunks.append(read_text_preview(p, max_lines=max_lines))
    return "".join(chunks)

def table_files_block(title, paths, max_rows=12):
    if not paths:
        return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
    chunks = [f"<h3>{html.escape(title)}</h3>"]
    for p in paths:
        chunks.append(f"<h4>{html.escape(p)}</h4>")
        chunks.append(read_table_preview(p, max_rows=max_rows))
    return "".join(chunks)

def image_files_block(title, paths):
    if not paths:
        return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
    chunks = [f"<h3>{html.escape(title)}</h3>"]
    for p in paths:
        asset = stage_asset(p)
        if asset is None:
            continue
        chunks.append(f"<h4>{html.escape(p)}</h4>")
        chunks.append(f'<img src="{html.escape(asset)}" alt="{html.escape(p)}" style="max-width: 1200px; width: 100%; border: 1px solid #ddd; margin-bottom: 14px;" />')
    return "".join(chunks)

def image_gallery_block(title, paths):
    if not paths:
        return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
    chunks = [f"<h3>{html.escape(title)}</h3>", '<div class="img-grid">']
    for p in paths:
        asset = stage_asset(p)
        if asset is None:
            continue
        chunks.append(
            '<figure class="img-card">'
            f'<img src="{html.escape(asset)}" alt="{html.escape(p)}" />'
            f'<figcaption>{html.escape(p)}</figcaption>'
            '</figure>'
        )
    chunks.append("</div>")
    return "".join(chunks)

def sample_from_report_path(rel_path):
    bits = rel_path.split("/")
    if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
        return bits[1]
    if len(bits) >= 3 and bits[0] == "ATAC":
        return bits[1]
    return None

def path_matches_sample(rel_path, sample):
    if not rel_path or not sample:
        return False
    bits = rel_path.split("/")
    if sample in bits:
        return True
    base = os.path.basename(rel_path)
    patterns = [
        f"{sample}.",
        f"{sample}_",
        f"{sample}-",
        f"{sample}R1",
        f"{sample}R2",
    ]
    return any(tok in base for tok in patterns)

def parse_starsolo_log_metrics(rel_path):
    wanted = [
        "Number of input reads",
        "Uniquely mapped reads %",
        "% of reads mapped to multiple loci",
        "% of reads unmapped: other",
        "% of reads unmapped: too short",
    ]
    out = {k: "" for k in wanted}
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return out
    try:
        with open(abs_path, "r", errors="replace") as fh:
            for line in fh:
                if "|" not in line:
                    continue
                left, right = line.split("|", 1)
                key = left.strip()
                val = right.strip()
                if key in out:
                    out[key] = val
    except Exception:
        return out
    return out

def starsolo_summary_table(log_paths):
    if not log_paths:
        return "<p><em>No STARsolo Log.final.out files found.</em></p>"
    headers = [
        "Sample",
        "Number of input reads",
        "Uniquely mapped reads %",
        "% of reads mapped to multiple loci",
        "% of reads unmapped: other",
        "% of reads unmapped: too short",
    ]
    rows = []
    for p in sorted(log_paths):
        bits = p.split("/")
        if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
            sample = bits[1]
        else:
            sample = os.path.basename(os.path.dirname(p))
        m = parse_starsolo_log_metrics(p)
        rows.append([
            sample,
            m["Number of input reads"],
            m["Uniquely mapped reads %"],
            m["% of reads mapped to multiple loci"],
            m["% of reads unmapped: other"],
            m["% of reads unmapped: too short"],
        ])
    cells = []
    cells.append("<tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + "</tr>")
    for row in rows:
        cells.append("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in row) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def _parse_pct(value):
    if value is None:
        return 0.0
    s = str(value).strip().replace("%", "").replace(",", "")
    try:
        return float(s)
    except Exception:
        return 0.0

def _parse_number(value):
    if value is None:
        return None
    s = str(value).strip().replace(",", "").replace("%", "")
    if s == "":
        return None
    try:
        return float(s)
    except Exception:
        return None

def _fmt_int(v):
    try:
        return f"{int(round(v)):,}"
    except Exception:
        return "N/A"

def _fmt_float(v, nd=1, suffix=""):
    try:
        return f"{float(v):.{nd}f}{suffix}"
    except Exception:
        return "N/A"

def load_experimental_groups(sample_barcode_path):
    groups = {}
    if not sample_barcode_path or not os.path.isfile(sample_barcode_path):
        return groups
    try:
        with open(sample_barcode_path, "r", errors="replace") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                cols = [c.strip() for c in re.split(r"\t|,", line)]
                if len(cols) < 4:
                    continue
                sample = cols[0]
                group = cols[3]
                if not sample or not group:
                    continue
                sample_l = sample.lower()
                group_l = group.lower()
                if sample_l in ("sample", "sample_name") and group_l in ("group", "experimental_group", "condition"):
                    continue
                groups[sample] = group
    except Exception:
        return {}
    return groups

def build_atac_cells_by_sample(atac_cell_summary_paths):
    atac_cells_by_sample = {}
    for p in sorted(atac_cell_summary_paths):
        s = atac_sample_key(p)
        if not s:
            continue
        atac_cells_by_sample[s] = parse_atac_cell_summary(p)
    return atac_cells_by_sample

def estimated_cells_for_sample(sample, summary_by_sample, atac_cells_by_sample):
    sy = summary_by_sample.get(sample, {})
    enc = sy.get("Estimated Number of Cells", "")
    atac_cells = atac_cells_by_sample.get(sample, {})
    if atac_cells:
        enc = (
            atac_cells.get("EstimatedCellsPreDedup", "")
            or atac_cells.get("EstimatedCells", "")
            or enc
        )
    return enc

def starsolo_metrics_by_sample(log_paths):
    out = {}
    for p in sorted(log_paths):
        sample = sample_from_report_path(p) or os.path.basename(os.path.dirname(p))
        out[sample] = parse_starsolo_log_metrics(p)
    return out

def summary_metrics_by_sample(paths):
    out = {}
    for p in sorted(paths):
        sample = sample_from_report_path(p) or os.path.basename(os.path.dirname(p))
        out[sample] = parse_summary_csv_metrics(p)
    return out

def overview_cards_html(sample_names, starsolo_by_sample, summary_by_sample):
    est_cells_vals = []
    for s in sample_names:
        sy = summary_by_sample.get(s, {})
        est = _parse_number(sy.get("Estimated Number of Cells", ""))
        if est is not None:
            est_cells_vals.append(est)
    n_samples = len(sample_names)
    total_est_cells = sum(est_cells_vals) if est_cells_vals else None
    generated = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    cards = []
    cards.append('<div class="kpi-grid">')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Samples</div><div class="kpi-value">{_fmt_int(n_samples)}</div></div>')
    cards.append(f'<div class="kpi-card"><div class="kpi-label">Total Estimated Cells - RNA</div><div class="kpi-value">{_fmt_int(total_est_cells)}</div></div>')
    cards.append('</div>')
    cards.append(f'<p class="meta-line"><strong>Generated:</strong> {html.escape(generated)} | <strong>Project:</strong> <code>{html.escape(proj)}</code></p>')
    return "".join(cards)

def sample_directory_table(sample_names, starsolo_by_sample, summary_by_sample, atac_cell_summary_paths):
    if not sample_names:
        return "<p><em>No samples detected.</em></p>"
    atac_cells_by_sample = build_atac_cells_by_sample(atac_cell_summary_paths)
    rows = []
    rows.append("<tr><th>Sample</th><th>Estimated Number of Cells</th><th>Sample report</th></tr>")
    for s in sample_names:
        enc = estimated_cells_for_sample(s, summary_by_sample, atac_cells_by_sample)
        link = f'QC_Report/{s}/index.html'
        rows.append(
            "<tr>"
            f"<td>{html.escape(s)}</td>"
            f"<td>{html.escape(enc)}</td>"
            f'<td><a href="{html.escape(link)}">Open</a></td>'
            "</tr>"
        )
    return "<table>" + "".join(rows) + "</table>"

def sample_combined_cell_cards(sample_names, summary_by_sample, atac_cell_summary_paths, sample_to_group):
    if not sample_names:
        return "<p><em>No samples detected.</em></p>"
    atac_cells_by_sample = build_atac_cells_by_sample(atac_cell_summary_paths)
    sample_estimates = {}
    for s in sample_names:
        sample_estimates[s] = _parse_number(estimated_cells_for_sample(s, summary_by_sample, atac_cells_by_sample))

    group_totals = {}
    if sample_to_group:
        members_by_group = {}
        for s in sample_names:
            g = sample_to_group.get(s, "")
            if g:
                members_by_group.setdefault(g, []).append(s)
        for g, members in members_by_group.items():
            vals = [sample_estimates.get(m) for m in members if sample_estimates.get(m) is not None]
            if vals:
                group_totals[g] = sum(vals)

    cards = ['<div class="kpi-grid">']
    for s in sorted(sample_names):
        g = sample_to_group.get(s, "")
        val = sample_estimates.get(s)
        label = s
        if g and g in group_totals:
            val = group_totals[g]
            label = f"{s} ({g})"
        cards.append(
            '<div class="kpi-card">'
            f'<div class="kpi-label">{html.escape(label)}</div>'
            f'<div class="kpi-value">{html.escape(_fmt_int(val))}</div>'
            '</div>'
        )
    cards.append('</div>')
    return "".join(cards)

def sample_sidebar_links(sample_names):
    if not sample_names:
        return "<p><em>No sample reports</em></p>"
    chunks = ['<div class="sample-side-title">Sample Reports</div>', '<div class="sample-side-links">']
    for s in sample_names:
        link = f'./QC_Report/{s}/index.html'
        chunks.append(f'<a href="{html.escape(link)}">{html.escape(s)}</a>')
    chunks.append('</div>')
    return "".join(chunks)

def alignment_summary_chart(log_paths):
    if not log_paths:
        return "<p><em>No STARsolo Log.final.out files found.</em></p>"
    rows = []
    for p in sorted(log_paths):
        bits = p.split("/")
        if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
            sample = bits[1]
        else:
            sample = os.path.basename(os.path.dirname(p))
        m = parse_starsolo_log_metrics(p)
        unique = _parse_pct(m.get("Uniquely mapped reads %", ""))
        multi = _parse_pct(m.get("% of reads mapped to multiple loci", ""))
        unmapped = max(0.0, 100.0 - unique - multi)
        total = unique + multi + unmapped
        if total <= 0:
            unique_w = multi_w = unmapped_w = 0.0
        else:
            unique_w = (unique / total) * 100.0
            multi_w = (multi / total) * 100.0
            unmapped_w = (unmapped / total) * 100.0
        rows.append((sample, unique, multi, unmapped, unique_w, multi_w, unmapped_w))

    chunks = []
    chunks.append('<div class="align-legend">')
    chunks.append('<span><i class="swatch unique"></i>Uniquely Mapped</span>')
    chunks.append('<span><i class="swatch multi"></i>Multi-mapped</span>')
    chunks.append('<span><i class="swatch unmapped"></i>Unmapped</span>')
    chunks.append('</div>')
    chunks.append('<div class="align-chart">')
    for sample, unique, multi, unmapped, unique_w, multi_w, unmapped_w in rows:
        chunks.append('<div class="align-row">')
        chunks.append(f'<div class="align-sample">{html.escape(sample)}</div>')
        chunks.append('<div class="align-bar">')
        chunks.append(f'<span class="seg unique" style="width:{unique_w:.4f}%"></span>')
        chunks.append(f'<span class="seg multi" style="width:{multi_w:.4f}%"></span>')
        chunks.append(f'<span class="seg unmapped" style="width:{unmapped_w:.4f}%"></span>')
        chunks.append('</div>')
        chunks.append(
            f'<div class="align-values">U: {unique:.2f}% | M: {multi:.2f}% | Un: {unmapped:.2f}%</div>'
        )
        chunks.append('</div>')
    chunks.append('</div>')
    return "".join(chunks)

def parse_barcodes_stats_metrics(rel_path):
    metrics = {}
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return metrics
    try:
        with open(abs_path, "r", errors="replace") as fh:
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                key = None
                val = None
                if "\t" in line:
                    parts = [x.strip() for x in line.split("\t") if x.strip()]
                    if len(parts) >= 2:
                        key = parts[0]
                        val = parts[-1]
                elif "|" in line:
                    left, right = line.split("|", 1)
                    key = left.strip()
                    val = right.strip()
                elif ":" in line:
                    left, right = line.split(":", 1)
                    key = left.strip()
                    val = right.strip()
                else:
                    parts = re.split(r"\s{2,}", line)
                    if len(parts) >= 2:
                        key = parts[0].strip()
                        val = parts[-1].strip()
                if key and val:
                    metrics[key] = val
    except Exception:
        return {}
    return metrics

def barcodes_stats_summary_table(paths):
    if not paths:
        return "<p><em>No Barcodes.stats files found.</em></p>"
    sample_metrics = []
    metric_order = []
    metric_seen = set()
    for p in sorted(paths):
        bits = p.split("/")
        if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
            sample = bits[1]
        else:
            sample = os.path.basename(os.path.dirname(p))
        metrics = parse_barcodes_stats_metrics(p)
        sample_metrics.append((sample, metrics))
        for k in metrics.keys():
            if k not in metric_seen:
                metric_seen.add(k)
                metric_order.append(k)
    sample_names = [s for s, _ in sample_metrics]
    headers = ["Metric"] + sample_names
    cells = []
    cells.append("<tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + "</tr>")
    for metric in metric_order:
        row = [metric] + [metrics.get(metric, "") for _, metrics in sample_metrics]
        cells.append("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in row) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def parse_summary_csv_metrics(rel_path):
    metrics = {}
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return metrics
    try:
        with open(abs_path, newline="") as fh:
            reader = csv.reader(fh)
            for row in reader:
                if not row:
                    continue
                key = row[0].strip()
                val = (row[1].strip() if len(row) > 1 else "")
                if key:
                    metrics[key] = val
    except Exception:
        return {}
    return metrics

def summary_csv_key_metrics_table(paths):
    wanted_metrics = [
        "Number of Reads",
        "Reads With Valid Barcodes",
        "Reads Mapped to Genome: Unique",
        "Estimated Number of Cells",
        "Median Reads per Cell",
        "Median UMI per Cell",
        "Total GeneFull Detected",
    ]
    if not paths:
        return "<p><em>No Summary.csv files found.</em></p>"
    sample_metrics = []
    for p in sorted(paths):
        bits = p.split("/")
        if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
            sample = bits[1]
        else:
            sample = os.path.basename(os.path.dirname(p))
        metrics = parse_summary_csv_metrics(p)
        sample_metrics.append((sample, metrics))
    headers = ["Metric"] + [s for s, _ in sample_metrics]
    cells = []
    cells.append("<tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + "</tr>")
    for metric in wanted_metrics:
        row_html = [f"<td>{html.escape(metric)}</td>"]
        for _, m in sample_metrics:
            val = m.get(metric, "")
            if metric == "Reads Mapped to Genome: Unique":
                pct = _parse_pct(val)
                if pct >= 85.0:
                    cls = "qc-good"
                elif pct >= 70.0:
                    cls = "qc-warn"
                else:
                    cls = "qc-bad"
                row_html.append(f'<td class="{cls}">{html.escape(str(val))}</td>')
            else:
                row_html.append(f"<td>{html.escape(str(val))}</td>")
        cells.append("<tr>" + "".join(row_html) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def parse_flagstat_metrics(rel_path):
    out = {
        "in_total": "",
        "mapped_pct": "",
        "properly_paired_pct": "",
        "duplicates_pct": "",
        "singletons_pct": "",
        "mate_diff_chr_count": "",
    }
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return out
    try:
        with open(abs_path, "r", errors="replace") as fh:
            for line in fh:
                s = line.strip()
                if " in total" in s and "+" in s:
                    out["in_total"] = s.split("+", 1)[0].strip()
                elif " mapped (" in s:
                    m = re.search(r'\(([^)]*%)', s)
                    out["mapped_pct"] = m.group(1) if m else ""
                elif " properly paired (" in s:
                    m = re.search(r'\(([^)]*%)', s)
                    out["properly_paired_pct"] = m.group(1) if m else ""
                elif " duplicates" in s and "(" in s:
                    m = re.search(r'\(([^)]*%)', s)
                    out["duplicates_pct"] = m.group(1) if m else ""
                elif " singletons (" in s:
                    m = re.search(r'\(([^)]*%)', s)
                    out["singletons_pct"] = m.group(1) if m else ""
                elif " with mate mapped to a different chr" in s and "(" not in s and "+" in s:
                    out["mate_diff_chr_count"] = s.split("+", 1)[0].strip()
    except Exception:
        return out
    return out

def parse_idxstats_metrics(rel_path):
    out = {
        "total_mapped": 0,
        "mt_mapped": 0,
        "autosomal_mapped": 0,
        "x_mapped": 0,
        "y_mapped": 0,
    }
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return out
    try:
        with open(abs_path, "r", errors="replace") as fh:
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 3:
                    continue
                chrom = cols[0].strip()
                try:
                    mapped = int(cols[2])
                except Exception:
                    continue
                out["total_mapped"] += mapped
                c = chrom.upper()
                if c in ("MT", "CHRM"):
                    out["mt_mapped"] += mapped
                elif c == "X" or c == "CHRX":
                    out["x_mapped"] += mapped
                elif c == "Y" or c == "CHRY":
                    out["y_mapped"] += mapped
                elif re.fullmatch(r'(CHR)?([1-9]|1[0-9]|2[0-2])', c):
                    out["autosomal_mapped"] += mapped
    except Exception:
        return out
    return out

def parse_atac_cell_summary(rel_path):
    out = {}
    abs_path = os.path.join(proj, rel_path)
    if not os.path.isfile(abs_path):
        return out
    try:
        with open(abs_path, "r", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    continue
                cols = line.rstrip("\n").split("\t", 1)
                if len(cols) == 2:
                    out[cols[0].strip()] = cols[1].strip()
    except Exception:
        return out
    return out

def parse_qc_metric_tsv(rel_path):
    return parse_atac_cell_summary(rel_path)

def atac_sample_key(rel_path):
    s = sample_from_report_path(rel_path)
    if s:
        return s
    b = os.path.basename(rel_path)
    m = re.match(r"(.+?)\.q30(?:\.(?:rmdup|mapped|possort))?\.(?:flagstat|idxstats|stats)\.txt\$", b)
    if m:
        return m.group(1)
    m2 = re.match(r"(.+?)\.atac_cells\.(?:summary|counts)\.tsv\$", b)
    if m2:
        return m2.group(1)
    return os.path.basename(os.path.dirname(rel_path))

def sample_name_matches(sample, candidate):
    if not sample or not candidate:
        return False
    s = str(sample).strip()
    c = str(candidate).strip()
    if c == s:
        return True
    s_l = s.lower()
    c_l = c.lower()
    if c_l == s_l:
        return True

    # MultiQC often stores full paths / filenames; compare against basename too.
    c_base = os.path.basename(c_l)

    def _norm(x):
        y = x
        y = y.replace(".fastq.gz", "").replace(".fq.gz", "")
        y = y.replace(".bam", "").replace(".txt", "")
        y = y.replace(".q30.rmdup.sorted", "").replace(".q30.rmdup", "")
        y = y.replace(".q30.possort", "").replace(".q30.mapped", "")
        y = y.replace(".flagstat", "").replace(".idxstats", "").replace(".stats", "")
        return y

    s_n = _norm(s_l)
    c_n = _norm(c_l)
    c_base_n = _norm(c_base)

    # Allow exact and containment matches after normalization.
    if s_n == c_n or s_n == c_base_n:
        return True
    if s_n in c_n or s_n in c_base_n:
        return True
    if c_n in s_n or c_base_n in s_n:
        return True
    return False

def sample_matches_atac_path(sample, rel_path):
    key = atac_sample_key(rel_path)
    return sample_name_matches(sample, key)

def multiqc_insert_size_for_sample(rel_paths, sample):
    pre_insert = ""
    post_insert = ""
    for rel in rel_paths:
        abs_path = os.path.join(proj, rel)
        if not os.path.isfile(abs_path):
            continue
        try:
            with open(abs_path, newline="") as fh:
                reader = csv.reader(fh, delimiter="\t")
                all_rows = list(reader)
        except Exception:
            continue
        if len(all_rows) < 2:
            continue
        header = all_rows[0]
        try:
            ins_idx = header.index("samtools_stats-insert_size_average")
        except ValueError:
            continue
        for r in all_rows[1:]:
            if not r:
                continue
            name = r[0] if len(r) > 0 else ""
            if sample_name_matches(sample, name):
                ins = (r[ins_idx] if ins_idx < len(r) else "").strip()
                if not ins:
                    continue
                name_l = str(name).lower()
                if ".q30.rmdup" in name_l:
                    post_insert = ins
                elif ".q30" in name_l:
                    pre_insert = ins
    return [pre_insert, post_insert]

def _fmt_pct(v):
    try:
        return f"{float(v):.2f}%"
    except Exception:
        return ""

def _to_int(s):
    try:
        return int(str(s).replace(",", "").strip())
    except Exception:
        return None

def atac_alignment_summary_table(flagstat_paths, idxstats_paths):
    if not flagstat_paths and not idxstats_paths:
        return "<p><em>No ATAC alignment metrics found.</em></p>"
    flag_by_sample = { (sample_from_report_path(p) or os.path.basename(os.path.dirname(p))): p for p in flagstat_paths }
    idx_by_sample = { (sample_from_report_path(p) or os.path.basename(os.path.dirname(p))): p for p in idxstats_paths }
    samples = sorted(set(list(flag_by_sample.keys()) + list(idx_by_sample.keys())))
    rows = []
    for sample in samples:
        f = parse_flagstat_metrics(flag_by_sample[sample]) if sample in flag_by_sample else {}
        i = parse_idxstats_metrics(idx_by_sample[sample]) if sample in idx_by_sample else {}
        total_mapped = i.get("total_mapped", 0)
        mt_frac = (100.0 * i.get("mt_mapped", 0) / total_mapped) if total_mapped > 0 else 0.0
        rows.append([
            sample,
            f.get("in_total", ""),
            f.get("properly_paired_pct", ""),
            f.get("singletons_pct", ""),
            f.get("mate_diff_chr_count", ""),
            _fmt_pct(mt_frac),
            f"{i.get('autosomal_mapped', 0):,}",
            f"{i.get('x_mapped', 0):,}",
            f"{i.get('y_mapped', 0):,}",
        ])
    cells = []
    cells.append(
        "<tr>"
        "<th>Sample</th>"
        "<th>Total reads (flagstat)</th>"
        "<th>Properly paired %</th>"
        "<th>Singleton %</th>"
        "<th>Mate on different chr</th>"
        "<th>Mitochondrial fraction</th>"
        "<th>Autosomal mapped</th>"
        "<th>X mapped</th>"
        "<th>Y mapped</th>"
        "</tr>"
    )
    for row in rows:
        cells.append("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in row) + "</tr>")
    return "<table>" + "".join(cells) + "</table>"

def atac_key_summary_table(flagstat_prededup_paths, idxstats_prededup_paths, flagstat_rmdup_paths, idxstats_rmdup_paths, atac_cell_summary_paths, atac_cbtag_qc_paths):
    samples = sorted(set(
        [atac_sample_key(p) for p in (flagstat_prededup_paths + idxstats_prededup_paths + flagstat_rmdup_paths + idxstats_rmdup_paths)]
    ))
    if not samples:
        return "<p><em>No ATAC alignment metrics found.</em></p>"

    pre_flag = { atac_sample_key(p): p for p in flagstat_prededup_paths }
    pre_idx = { atac_sample_key(p): p for p in idxstats_prededup_paths }
    post_flag = { atac_sample_key(p): p for p in flagstat_rmdup_paths }
    post_idx = { atac_sample_key(p): p for p in idxstats_rmdup_paths }
    cell_sum = { atac_sample_key(p): p for p in atac_cell_summary_paths }
    cbtag_qc = { atac_sample_key(p): p for p in atac_cbtag_qc_paths }

    by_sample = {}
    for s in samples:
        pf = parse_flagstat_metrics(pre_flag[s]) if s in pre_flag else {}
        pi = parse_idxstats_metrics(pre_idx[s]) if s in pre_idx else {}
        rf = parse_flagstat_metrics(post_flag[s]) if s in post_flag else {}
        ri = parse_idxstats_metrics(post_idx[s]) if s in post_idx else {}
        cs = parse_atac_cell_summary(cell_sum[s]) if s in cell_sum else {}
        cq = parse_qc_metric_tsv(cbtag_qc[s]) if s in cbtag_qc else {}

        pre_total = _to_int(pf.get("in_total", ""))
        post_total = _to_int(rf.get("in_total", ""))
        retained_pct = (100.0 * post_total / pre_total) if (pre_total and post_total is not None and pre_total > 0) else None

        pre_mt = (100.0 * pi.get("mt_mapped", 0) / pi.get("total_mapped", 1)) if pi.get("total_mapped", 0) > 0 else 0.0
        post_mt = (100.0 * ri.get("mt_mapped", 0) / ri.get("total_mapped", 1)) if ri.get("total_mapped", 0) > 0 else 0.0

        by_sample[s] = {
            "EstimatedCellsPreDedup": cs.get("EstimatedCellsPreDedup", ""),
            "EstimatedCells": cs.get("EstimatedCells", ""),
            "MedianFragmentsPreDedup": cs.get("MedianFragmentsPerEstimatedCellPreDedup", ""),
            "MedianFragmentsPostDedup": cs.get("MedianFragmentsPerEstimatedCell", ""),
            "MedianTSSPreDedup": cs.get("MedianTSSEnrichmentPerEstimatedCellPreDedup", ""),
            "MedianTSSPostDedup": cs.get("MedianTSSEnrichmentPerEstimatedCell", ""),
            "MedianReadsInTSSPreDedup": cs.get("MedianReadsInTSSPerEstimatedCellPreDedup", ""),
            "MedianReadsInTSSPostDedup": cs.get("MedianReadsInTSSPerEstimatedCell", ""),
            "MedianReadsInPromoterPreDedup": cs.get("MedianReadsInPromoterPerEstimatedCellPreDedup", ""),
            "MedianReadsInPromoterPostDedup": cs.get("MedianReadsInPromoterPerEstimatedCell", ""),
            "MedianPromoterRatioPreDedup": cs.get("MedianPromoterRatioPerEstimatedCellPreDedup", ""),
            "MedianPromoterRatioPostDedup": cs.get("MedianPromoterRatioPerEstimatedCell", ""),
            "MedianReadsInBlacklistPreDedup": cs.get("MedianReadsInBlacklistPerEstimatedCellPreDedup", ""),
            "MedianReadsInBlacklistPostDedup": cs.get("MedianReadsInBlacklistPerEstimatedCell", ""),
            "MedianBlacklistRatioPreDedup": cs.get("MedianBlacklistRatioPerEstimatedCellPreDedup", ""),
            "MedianBlacklistRatioPostDedup": cs.get("MedianBlacklistRatioPerEstimatedCell", ""),
            "CBTaggedPct": cq.get("CBTaggedPct", ""),
            "ReadsPreDedup": pf.get("in_total", ""),
            "MTPctPreDedup": _fmt_pct(pre_mt),
            "ReadsPostDedup": rf.get("in_total", ""),
            "MTPctPostDedup": _fmt_pct(post_mt),
            "ReadsRetainedPct": _fmt_pct(retained_pct) if retained_pct is not None else "",
        }

    row_defs = [
        ("Estimated cells (ArchR, pre-dedup)", "EstimatedCellsPreDedup"),
        ("Estimated cells (ArchR, post-dedup)", "EstimatedCells"),
        ("Total reads (flagstat, pre-dedup)", "ReadsPreDedup"),
        ("Total reads (flagstat, post-dedup)", "ReadsPostDedup"),
        ("Reads retained after dedup (%)", "ReadsRetainedPct"),
        ("Median nFrags per cell (ArchR, pre-dedup)", "MedianFragmentsPreDedup"),
        ("Median nFrags per cell (ArchR, post-dedup)", "MedianFragmentsPostDedup"),
        ("Median TSS enrichment (ArchR, pre-dedup)", "MedianTSSPreDedup"),
        ("Median TSS enrichment (ArchR, post-dedup)", "MedianTSSPostDedup"),
    ]

    cells = []
    cells.append('<tr><th class="atac-pivot-metric">Metric</th>')
    for s in samples:
        cells.append(f'<th>{html.escape(str(s))}</th>')
    cells.append("</tr>")
    for label, key in row_defs:
        cells.append(f'<tr><td class="atac-pivot-metric">{html.escape(label)}</td>')
        for s in samples:
            v = by_sample.get(s, {}).get(key, "")
            cells.append(f"<td>{html.escape(str(v))}</td>")
        cells.append("</tr>")
    return '<div class="atac-pivot-wrap"><table class="atac-key-pivot">' + "".join(cells) + "</table></div>"

demux_total = rel_list("demux/*.total_number_reads.tsv")
demux_stats = sorted(set(
    rel_list("demux/SHARE-seq.demultiplex.stats.tsv")
    + rel_list_recursive("demux/**/SHARE-seq.demultiplex.stats.tsv")
))
fastqc_html = sorted(set(
    rel_list("fastqc_demux/*/*_fastqc.html") +
    rel_list("fastqc_demux/*/fastqc_demux/*_fastqc.html") +
    rel_list("fastqc_demux/*_fastqc.html")
))
fastqc_trimmed_html = sorted(set(
    rel_list("fastqc_trimmed/*_fastqc.html") +
    rel_list_recursive("fastqc_trimmed/**/*_fastqc.html")
))

starsolo_logs = sorted(set(
    p for root in active_starsolo_roots
    for p in rel_list(f"{root}/*/Log.final.out")
))
knee_plots = sorted(set(
    p for root in active_starsolo_roots
    for p in (rel_list(f"{root}/*/*_knee_plot.png") + rel_list_recursive(f"{root}/**/*_knee_plot.png"))
))
barcodes_stats = sorted(set(
    p for root in active_starsolo_roots
    for p in rel_list(f"{root}/*/Solo.out/Barcodes.stats")
))
summary_csv = sorted(set(
    p for root in active_starsolo_roots
    for p in rel_list(f"{root}/*/Solo.out/GeneFull/Summary.csv")
))
barnyard = sorted(set(
    p for root in active_starsolo_roots
    for p in rel_list(f"{root}/*/*collision_plot.png")
))
atac_flagstat_rmdup = rel_list("ATAC/*/*.q30.rmdup.flagstat.txt")
atac_idxstats_rmdup = rel_list("ATAC/*/*.q30.rmdup.idxstats.txt")
atac_stats_rmdup = rel_list("ATAC/*/*.q30.rmdup.stats.txt")
atac_flagstat_prededup = rel_list("ATAC/*/*.q30.mapped.flagstat.txt")
atac_idxstats_prededup = rel_list("ATAC/*/*.q30.mapped.idxstats.txt")
atac_stats_prededup = rel_list("ATAC/*/*.q30.mapped.stats.txt")
atac_cbtag_qc = rel_list("ATAC/*/*.cbtag_qc.tsv")
atac_cell_summary = sorted(set(
    rel_list("ATAC/*/*.atac_cells.summary.tsv") +
    rel_list_recursive("ATAC/**/*.atac_cells.summary.tsv")
))
atac_multiqc_report = sorted(set(
    rel_list("multiqc_atac/ATAC_MultiQC.html") +
    rel_list("ATAC_MultiQC.html")
))
atac_multiqc_tables = sorted(set(
    rel_list("multiqc_atac/ATAC_MultiQC_data/multiqc_general_stats.txt") +
    rel_list("multiqc_atac/ATAC_MultiQC_data/multiqc_samtools_flagstat.txt") +
    rel_list("multiqc_atac/ATAC_MultiQC_data/multiqc_samtools_stats.txt") +
    rel_list("ATAC_MultiQC_data/multiqc_general_stats.txt") +
    rel_list("ATAC_MultiQC_data/multiqc_samtools_flagstat.txt") +
    rel_list("ATAC_MultiQC_data/multiqc_samtools_stats.txt")
))
sample_names = sorted(set(
    [sample_from_report_path(p) for p in (starsolo_logs + barcodes_stats + summary_csv + knee_plots + barnyard + atac_flagstat_rmdup + atac_idxstats_rmdup + atac_stats_rmdup + atac_flagstat_prededup + atac_idxstats_prededup + atac_stats_prededup + atac_cbtag_qc) if sample_from_report_path(p)]
    + list(load_demux_sample_names(demux_stats))
))
starsolo_by_sample = starsolo_metrics_by_sample(starsolo_logs)
summary_by_sample = summary_metrics_by_sample(summary_csv)
sample_to_group = load_experimental_groups(sample_barcode_sheet)

parts = []
parts.append('''<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SHARE-seq QC Report</title>
  <style>
    body {
      font-family: "Inter", "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      font-size: 14px;
      color: #000000;
      margin: 24px auto;
      max-width: 1400px;
      line-height: 1.55;
      padding: 0 12px;
      background: #ffffff;
      scroll-padding-top: 86px;
    }
    h1, h2, h3, h4 {
      color: #000000;
      margin-top: 1.15em;
      margin-bottom: 0.45em;
      font-weight: 700;
      letter-spacing: 0.1px;
    }
    h1 { font-size: 30px; margin-top: 0.3em; }
    h2 { font-size: 21px; border-bottom: 1px solid #dfe1df; padding-bottom: 4px; }
    h3 { font-size: 16px; font-weight: 700; }
    h4 { font-size: 13px; font-weight: 600; color: #000000; }
    p { margin: 0.45em 0 0.8em 0; }
    em { color: #000000; }
    code { background: #ffffff; color: #000000; border: 1px solid #dfe1df; padding: 2px 5px; border-radius: 4px; font-size: 0.92em; }
    table {
      border-collapse: collapse;
      margin: 10px 0 20px 0;
      width: 100%;
      background: #fff;
      box-shadow: 0 1px 2px rgba(141, 0, 52, 0.08);
      display: table;
    }
    th, td {
      border: 1px solid #dfe1df;
      padding: 6px 10px;
      font-size: 12px;
      vertical-align: top;
      white-space: nowrap;
    }
    th {
      background: #ffffff;
      font-weight: 700;
      text-align: left;
      color: #000000;
    }
    td:first-child { font-weight: 600; color: #000000; }
    tr:nth-child(even) td { background: #ffffff; }
    .atac-pivot-wrap { overflow-x: auto; margin: 8px 0 14px 0; max-width: 100%; -webkit-overflow-scrolling: touch; }
    table.atac-key-pivot { width: max-content; max-width: none; }
    .atac-key-pivot th.atac-pivot-metric,
    .atac-key-pivot td.atac-pivot-metric { white-space: normal; min-width: 200px; max-width: min(360px, 40vw); }
    .atac-key-pivot th:not(.atac-pivot-metric),
    .atac-key-pivot td:not(.atac-pivot-metric) { white-space: nowrap; }
    .tabs {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      gap: 8px;
      padding-left: 10px;
      padding-right: 8px;
    }
    .tabs a {
      text-decoration: none;
      border: 1px solid rgba(255, 255, 255, 0.6);
      background: rgba(255, 255, 255, 0.14);
      color: #ffffff;
      padding: 7px 12px;
      border-radius: 6px;
      font-size: 12px;
      font-weight: 600;
    }
    .tabs a:hover { background: rgba(255, 255, 255, 0.26); border-color: #ffffff; }
    .top-banner {
      background: #d11947;
      color: #ffffff;
      border-radius: 10px;
      padding: 12px 14px 12px 14px;
      margin: 0 0 14px 0;
      box-shadow: 0 2px 6px rgba(141, 0, 52, 0.22);
    }
    .top-banner h1, .top-banner h2, .top-banner h3, .top-banner h4,
    .top-banner p, .top-banner em, .top-banner code { color: #ffffff; }
    .top-banner code {
      background: rgba(255, 255, 255, 0.14);
      border: 1px solid rgba(255, 255, 255, 0.45);
    }
    section {
      padding: 10px 14px 6px 14px;
      margin: 8px 0 14px 0;
      background: #ffffff;
      color: #000000;
      border: 1px solid #dfe1df;
      border-radius: 10px;
      box-shadow: 0 1px 2px rgba(141, 0, 52, 0.06);
      scroll-margin-top: 90px;
    }
    section h2, section h3, section h4, section p, section em, section a { color: #000000; }
    section h2 { border-bottom: 1px solid #dfe1df; }
    section code {
      background: #ffffff;
      border: 1px solid #dfe1df;
      color: #000000;
    }
    section table, section th, section td {
      background: #ffffff;
      color: #000000;
    }
    .img-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
      gap: 12px;
      margin: 8px 0 18px 0;
    }
    .img-card {
      margin: 0;
      border: 1px solid #dfe1df;
      border-radius: 6px;
      padding: 8px;
      background: #fff;
    }
    .img-card img {
      width: 100%;
      height: auto;
      max-height: 260px;
      object-fit: contain;
      display: block;
    }
    .img-card figcaption {
      font-size: 11px;
      margin-top: 6px;
      word-break: break-word;
    }
    .align-legend { display: flex; gap: 14px; flex-wrap: wrap; margin: 6px 0 10px 0; font-size: 12px; }
    .align-legend .swatch {
      display: inline-block;
      width: 10px;
      height: 10px;
      border-radius: 2px;
      margin-right: 5px;
      vertical-align: middle;
    }
    .align-chart { margin: 6px 0 16px 0; }
    .align-row {
      display: grid;
      grid-template-columns: 180px minmax(280px, 1fr) 260px;
      gap: 10px;
      align-items: center;
      margin-bottom: 8px;
    }
    .align-sample { font-size: 12px; word-break: break-word; }
    .align-bar {
      display: flex;
      height: 14px;
      border: 1px solid #d0d0d0;
      border-radius: 3px;
      overflow: hidden;
      background: #f5f5f5;
    }
    .align-bar .seg { height: 100%; display: block; }
    .align-values { font-size: 11px; color: #000000; font-weight: 600; }
    .unique { background: #2e7d32; }
    .multi { background: #1976d2; }
    .unmapped { background: #d32f2f; }
    section td.qc-good { background: #e6f4ea !important; color: #000000 !important; font-weight: 700; }
    section td.qc-warn { background: #fff8e1 !important; color: #000000 !important; font-weight: 700; }
    section td.qc-bad { background: #f7b8b8 !important; color: #000000 !important; font-weight: 700; }
    .kpi-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
      gap: 10px;
      margin: 8px 0 10px 0;
    }
    .kpi-card {
      border: 1px solid #dfe1df;
      border-radius: 8px;
      background: #ffffff;
      padding: 10px 12px;
    }
    .kpi-label {
      font-size: 12px;
      font-weight: 600;
      color: #000000;
      margin-bottom: 2px;
    }
    .kpi-value {
      font-size: 24px;
      font-weight: 700;
      color: #000000;
      line-height: 1.2;
    }
    .meta-line { margin-top: 4px; color: #000000; }
    .main-layout {
      display: grid;
      grid-template-columns: 230px minmax(0, 1fr);
      gap: 14px;
      align-items: start;
    }
    .sample-sidebar {
      position: sticky;
      top: 74px;
      align-self: start;
      background: #d11947;
      border: 1px solid #8d0034;
      border-radius: 10px;
      padding: 10px;
      color: #ffffff;
      box-shadow: 0 2px 6px rgba(141, 0, 52, 0.22);
      max-height: calc(100vh - 74px);
      overflow-y: auto;
    }
    .tabs-wrap {
      position: sticky;
      top: 0;
      z-index: 1000;
      background: #d11947;
      border-radius: 0 0 8px 8px;
      padding: 8px 0 6px 0;
      margin-bottom: 14px;
      box-shadow: 0 2px 4px rgba(141, 0, 52, 0.18);
    }
    .sample-side-title {
      font-size: 12px;
      font-weight: 700;
      color: #ffffff;
      margin-bottom: 8px;
      text-transform: uppercase;
      letter-spacing: 0.3px;
    }
    .sample-side-links {
      display: flex;
      flex-direction: column;
      gap: 6px;
    }
    .sample-side-links a {
      text-decoration: none;
      font-size: 12px;
      color: #ffffff;
      border: 1px solid rgba(255, 255, 255, 0.6);
      border-radius: 6px;
      padding: 6px 8px;
      background: rgba(255, 255, 255, 0.14);
      font-weight: 700;
    }
    .sample-side-links a:hover { background: rgba(255, 255, 255, 0.26); border-color: #ffffff; }
    .main-content { min-width: 0; }
    .starsolo-block { margin: 0 0 18px 0; }
    @media (max-width: 1100px) {
      .main-layout { grid-template-columns: 1fr; }
      .sample-sidebar { position: static; max-height: none; }
    }
  </style>
</head>
<body>
<div class="top-banner">
<h1>SHARE-seq QC and Visualization Report</h1>
<p>Generated from pipeline outputs in <code>__PROJECT_DIR__</code>.</p>
</div>
<div class="tabs-wrap">
<div class="tabs">
  <a href="#sec-overview">Overview</a>
  <a href="#sec-demux">Demultiplexing</a>
  <a href="#sec-fastqc">FastQC (Demultiplexed)</a>
  <a href="#sec-atac">ATAC QC</a>
  <a href="#sec-starsolo">RNA QC</a>
__BARNYARD_TAB__
</div>
</div>
<div class="main-layout">
<aside class="sample-sidebar">
__SAMPLE_SIDEBAR__
</aside>
<div class="main-content">
'''.replace("__PROJECT_DIR__", html.escape(proj)))
parts[-1] = parts[-1].replace("__SAMPLE_SIDEBAR__", sample_sidebar_links(sample_names))
parts[-1] = parts[-1].replace(
    "__BARNYARD_TAB__",
    '  <a href="#sec-barnyard">Hybrid Barnyard</a>' if species_model == 'hybrid' else ""
)

parts.append('<section id="sec-overview">')
parts.append("<h2>Overview</h2>")
parts.append(overview_cards_html(sample_names, starsolo_by_sample, summary_by_sample))
parts.append("<h3>Combined Cell Estimate by Sample</h3>")
parts.append(sample_combined_cell_cards(sample_names, summary_by_sample, atac_cell_summary, sample_to_group))
parts.append("<h3>Sample Directory</h3>")
parts.append(sample_directory_table(sample_names, starsolo_by_sample, summary_by_sample, atac_cell_summary))
parts.append("</section>")

parts.append('<section id="sec-demux">')
parts.append("<h2>Demultiplexing</h2>")
if demux_stats:
    parts.append("<h3>Demultiplex Stats</h3>")
    for dsp in demux_stats:
        if len(demux_stats) > 1:
            parts.append(f"<h4>{html.escape(dsp)}</h4>")
        parts.append(read_table_preview(dsp, max_rows=None))
else:
    parts.append("<p><em>No demultiplex stats file found.</em></p>")
parts.append("</section>")

parts.append('<section id="sec-fastqc">')
parts.append("<h2>FastQC (Demultiplexed)</h2>")
parts.append(links_block("fastqc_demux/<sample>/*_fastqc.html", fastqc_html))
parts.append("</section>")

parts.append('<section id="sec-atac">')
parts.append("<h2>ATAC QC</h2>")
parts.append(
    "<p><strong>What these files mean:</strong> "
    "<code>.q30.mapped</code> metrics are captured after MAPQ>=30 filtering and before duplicate removal. "
    "<code>.q30.rmdup</code> metrics are captured after duplicate removal. "
    "Use the difference to estimate duplicate burden and retained unique signal.</p>"
)
parts.append('<div class="starsolo-block"><h3>ATAC Key Summary by Sample</h3>')
parts.append(atac_key_summary_table(atac_flagstat_prededup, atac_idxstats_prededup, atac_flagstat_rmdup, atac_idxstats_rmdup, atac_cell_summary, atac_cbtag_qc))
parts.append("</div>")
parts.append("</section>")

parts.append('<section id="sec-starsolo">')
parts.append("<h2>RNA QC</h2>")
parts.append('<div class="starsolo-block"><h3>Summary.csv Key Metrics by Sample</h3>')
parts.append(summary_csv_key_metrics_table(summary_csv))
parts.append("</div>")
parts.append('<div class="starsolo-block"><h3>Alignment Summary Bar Chart</h3>')
parts.append(alignment_summary_chart(starsolo_logs))
parts.append("</div>")
parts.append('<div class="starsolo-block">')
parts.append(image_gallery_block("STARsolo*/<sample>/*_knee_plot.png", knee_plots))
parts.append("</div>")
parts.append("</section>")

if species_model == 'hybrid':
    parts.append('<section id="sec-barnyard">')
    parts.append("<h2>Hybrid Barnyard Plots</h2>")
    parts.append(image_files_block("STARsolo*/<sample>/*collision_plot.png", barnyard))
    parts.append("</section>")

parts.append("</div></div>")
parts.append("</body></html>")

with open(out_path, "w") as out:
    out.write("\n".join(parts))

def sample_from_starsolo_path(rel_path):
    bits = rel_path.split("/")
    if len(bits) >= 3 and bits[0] in ("STARsolo", "STARsolo_paired"):
        return bits[1]
    return None

def sample_from_atac_path(rel_path):
    bits = rel_path.split("/")
    if len(bits) >= 3 and bits[0] == "ATAC":
        return bits[1]
    return None

all_sample_candidates = set()
for p in starsolo_logs + knee_plots + barcodes_stats + summary_csv + barnyard + atac_flagstat_rmdup + atac_idxstats_rmdup + atac_stats_rmdup + atac_flagstat_prededup + atac_idxstats_prededup + atac_stats_prededup + atac_cbtag_qc:
    s = sample_from_starsolo_path(p)
    if s:
        all_sample_candidates.add(s)
    a = sample_from_atac_path(p)
    if a:
        all_sample_candidates.add(a)
all_sample_candidates.update(load_demux_sample_names(demux_stats))

for sample in sorted(all_sample_candidates):
    sample_root = os.path.join(per_sample_dir, sample)
    sample_assets = os.path.join(sample_root, "assets")
    os.makedirs(sample_assets, exist_ok=True)

    sample_logs = [p for p in starsolo_logs if sample_from_starsolo_path(p) == sample]
    sample_knee = [p for p in knee_plots if sample_from_starsolo_path(p) == sample]
    sample_barcodes = [p for p in barcodes_stats if sample_from_starsolo_path(p) == sample]
    sample_summary = [p for p in summary_csv if sample_from_starsolo_path(p) == sample]
    sample_barnyard = [p for p in barnyard if sample_from_starsolo_path(p) == sample]
    sample_fastqc_demux = [p for p in fastqc_html if path_matches_sample(p, sample)]
    sample_fastqc_trimmed = [p for p in fastqc_trimmed_html if path_matches_sample(p, sample)]
    sample_atac_flagstat_rmdup = [p for p in atac_flagstat_rmdup if sample_matches_atac_path(sample, p)]
    sample_atac_idxstats_rmdup = [p for p in atac_idxstats_rmdup if sample_matches_atac_path(sample, p)]
    sample_atac_stats_rmdup = [p for p in atac_stats_rmdup if sample_matches_atac_path(sample, p)]
    sample_atac_flagstat_prededup = [p for p in atac_flagstat_prededup if sample_matches_atac_path(sample, p)]
    sample_atac_idxstats_prededup = [p for p in atac_idxstats_prededup if sample_matches_atac_path(sample, p)]
    sample_atac_stats_prededup = [p for p in atac_stats_prededup if sample_matches_atac_path(sample, p)]
    sample_atac_cbtag_qc = [p for p in atac_cbtag_qc if sample_matches_atac_path(sample, p)]
    sample_atac_cells = [p for p in atac_cell_summary if sample_matches_atac_path(sample, p)]
    sample_atac_multiqc_tables = atac_multiqc_tables
    has_atac = bool(sample_atac_flagstat_rmdup or sample_atac_idxstats_rmdup or sample_atac_stats_rmdup or sample_atac_flagstat_prededup or sample_atac_idxstats_prededup or sample_atac_stats_prededup or sample_atac_cbtag_qc)
    has_rna = bool(sample_logs or sample_knee or sample_barcodes or sample_summary or sample_barnyard)
    sample_tabs = [("sec-demux", "Demultiplexing")]
    if has_atac:
        sample_tabs.append(("sec-atac", "ATAC QC"))
    if has_rna:
        sample_tabs.append(("sec-rna", "RNA QC"))

    def sample_image_block(title, paths):
        if not paths:
            return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
        chunks = [f"<h3>{html.escape(title)}</h3>"]
        for p in paths:
            asset = stage_asset_in(p, sample_assets)
            if asset is None:
                continue
            rel_asset = os.path.relpath(asset, sample_root)
            chunks.append(f"<h4>{html.escape(p)}</h4>")
            chunks.append(f'<img src="{html.escape(rel_asset)}" alt="{html.escape(p)}" style="max-width: 1200px; width: 100%; border: 1px solid #ddd; margin-bottom: 14px;" />')
        return "".join(chunks)

    def sample_links_block(title, paths):
        if not paths:
            return f"<h3>{html.escape(title)}</h3><p><em>No files found.</em></p>"
        items = []
        for p in paths:
            asset = stage_asset_in(p, sample_assets)
            if asset is None:
                continue
            rel_asset = os.path.relpath(asset, sample_root)
            items.append(f'<li><a href="{html.escape(rel_asset)}">{html.escape(p)}</a></li>')
        if not items:
            return f"<h3>{html.escape(title)}</h3><p><em>No readable files found.</em></p>"
        return f"<h3>{html.escape(title)}</h3><ul>{''.join(items)}</ul>"

    def prefer_path(paths, must_contain):
        if not paths:
            return None
        preferred = [p for p in paths if must_contain in p]
        pool = preferred if preferred else paths
        return sorted(pool)[0] if pool else None

    def filter_paths(paths, must_contain):
        if not paths:
            return []
        return sorted([p for p in paths if must_contain in p])

    def atac_pre_post_table(sample):
        pre_flag_path = prefer_path(sample_atac_flagstat_prededup, ".q30.mapped.flagstat.txt")
        pre_idx_path = prefer_path(sample_atac_idxstats_prededup, ".q30.mapped.idxstats.txt")
        post_flag_path = prefer_path(sample_atac_flagstat_rmdup, ".q30.rmdup.flagstat.txt")
        post_idx_path = prefer_path(sample_atac_idxstats_rmdup, ".q30.rmdup.idxstats.txt")

        pre_flag = parse_flagstat_metrics(pre_flag_path) if pre_flag_path else {}
        pre_idx = parse_idxstats_metrics(pre_idx_path) if pre_idx_path else {}
        post_flag = parse_flagstat_metrics(post_flag_path) if post_flag_path else {}
        post_idx = parse_idxstats_metrics(post_idx_path) if post_idx_path else {}
        pre_insert_avg, post_insert_avg = multiqc_insert_size_for_sample(sample_atac_multiqc_tables, sample)
        cs = parse_atac_cell_summary(sample_atac_cells[0]) if sample_atac_cells else {}
        cq = parse_qc_metric_tsv(sample_atac_cbtag_qc[0]) if sample_atac_cbtag_qc else {}

        pre_total = _to_int(pre_flag.get("in_total", ""))
        post_total = _to_int(post_flag.get("in_total", ""))
        retained_pct = (100.0 * post_total / pre_total) if (pre_total and post_total is not None and pre_total > 0) else None
        removed_pct = (100.0 - retained_pct) if retained_pct is not None else None
        pre_mt = (100.0 * pre_idx.get("mt_mapped", 0) / pre_idx.get("total_mapped", 1)) if pre_idx.get("total_mapped", 0) > 0 else 0.0
        post_mt = (100.0 * post_idx.get("mt_mapped", 0) / post_idx.get("total_mapped", 1)) if post_idx.get("total_mapped", 0) > 0 else 0.0

        rows = [
            ("Estimated cells (ArchR)", cs.get("EstimatedCellsPreDedup", ""), cs.get("EstimatedCells", "")),
            ("Median nFrags (ArchR)", cs.get("MedianFragmentsPerEstimatedCellPreDedup", ""), cs.get("MedianFragmentsPerEstimatedCell", "")),
            ("Median TSS enrichment (ArchR)", cs.get("MedianTSSEnrichmentPerEstimatedCellPreDedup", ""), cs.get("MedianTSSEnrichmentPerEstimatedCell", "")),
            ("ReadsInTSS (ArchR)", cs.get("MedianReadsInTSSPerEstimatedCellPreDedup", ""), cs.get("MedianReadsInTSSPerEstimatedCell", "")),
            ("CB-tagged alignments %", cq.get("CBTaggedPct", ""), ""),
            ("CB-tagged alignments", cq.get("CBTaggedAlignments", ""), ""),
            ("Missing CB alignments", cq.get("MissingCBAlignments", ""), ""),
            ("ReadsInPromoter (ArchR)", cs.get("MedianReadsInPromoterPerEstimatedCellPreDedup", ""), cs.get("MedianReadsInPromoterPerEstimatedCell", "")),
            ("PromoterRatio (ArchR)", cs.get("MedianPromoterRatioPerEstimatedCellPreDedup", ""), cs.get("MedianPromoterRatioPerEstimatedCell", "")),
            ("ReadsInBlacklist (ArchR)", cs.get("MedianReadsInBlacklistPerEstimatedCellPreDedup", ""), cs.get("MedianReadsInBlacklistPerEstimatedCell", "")),
            ("BlacklistRatio (ArchR)", cs.get("MedianBlacklistRatioPerEstimatedCellPreDedup", ""), cs.get("MedianBlacklistRatioPerEstimatedCell", "")),
            ("Total reads", pre_flag.get("in_total", ""), post_flag.get("in_total", "")),
            ("Properly paired %", pre_flag.get("properly_paired_pct", ""), post_flag.get("properly_paired_pct", "")),
            ("Singleton %", pre_flag.get("singletons_pct", ""), post_flag.get("singletons_pct", "")),
            ("Mate on different chr", pre_flag.get("mate_diff_chr_count", ""), post_flag.get("mate_diff_chr_count", "")),
            ("Insert size avg", pre_insert_avg, post_insert_avg),
            ("Mitochondrial fraction", _fmt_pct(pre_mt), _fmt_pct(post_mt)),
            ("Autosomal mapped", f"{pre_idx.get('autosomal_mapped', 0):,}", f"{post_idx.get('autosomal_mapped', 0):,}"),
            ("X mapped", f"{pre_idx.get('x_mapped', 0):,}", f"{post_idx.get('x_mapped', 0):,}"),
            ("Y mapped", f"{pre_idx.get('y_mapped', 0):,}", f"{post_idx.get('y_mapped', 0):,}"),
            ("Reads retained", "", _fmt_pct(retained_pct) if retained_pct is not None else ""),
            ("Reads removed", "", _fmt_pct(removed_pct) if removed_pct is not None else ""),
        ]
        cells = ['<tr><th>Metric</th><th>Pre-dedup (.q30.mapped)</th><th>Post-dedup (.q30.rmdup)</th></tr>']
        for m, pre_v, post_v in rows:
            cells.append(f"<tr><td>{html.escape(str(m))}</td><td>{html.escape(str(pre_v))}</td><td>{html.escape(str(post_v))}</td></tr>")
        return "<table>" + "".join(cells) + "</table>"

    sample_parts = []
    sample_parts.append(f'''<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SHARE-seq QC Report - {html.escape(sample)}</title>
  <style>
    body {{
      font-family: "Inter", "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      font-size: 14px;
      color: #000000;
      margin: 24px auto;
      max-width: 1200px;
      line-height: 1.55;
      padding: 0 12px;
      background: #ffffff;
      scroll-padding-top: 78px;
    }}
    h1, h2, h3 {{
      color: #000000;
      margin-top: 1.15em;
      margin-bottom: 0.45em;
      font-weight: 700;
      letter-spacing: 0.1px;
    }}
    h1 {{ font-size: 28px; margin-top: 0.3em; }}
    h2 {{ font-size: 20px; border-bottom: 1px solid #dfe1df; padding-bottom: 4px; }}
    h3 {{ font-size: 16px; }}
    code {{ background: #ffffff; color: #000000; border: 1px solid #dfe1df; padding: 2px 5px; border-radius: 4px; font-size: 0.92em; }}
    pre {{
      background: #ffffff;
      border: 1px solid #dfe1df;
      padding: 10px;
      overflow-x: auto;
      line-height: 1.45;
      font-size: 12px;
    }}
    table {{
      border-collapse: collapse;
      margin: 10px 0 20px 0;
      width: 100%;
      background: #fff;
      box-shadow: 0 1px 2px rgba(141, 0, 52, 0.08);
    }}
    th, td {{
      border: 1px solid #dfe1df;
      padding: 6px 10px;
      font-size: 12px;
      vertical-align: top;
    }}
    th {{ background: #ffffff; font-weight: 700; text-align: left; color: #000000; }}
    tr:nth-child(even) td {{ background: #ffffff; }}
    .top-banner {{
      background: #d11947;
      color: #ffffff;
      border-radius: 10px;
      padding: 12px 14px 12px 14px;
      margin: 0 0 14px 0;
      box-shadow: 0 2px 6px rgba(141, 0, 52, 0.22);
    }}
    .top-banner h1, .top-banner p, .top-banner code {{ color: #ffffff; }}
    .top-banner code {{
      background: rgba(255, 255, 255, 0.14);
      border: 1px solid rgba(255, 255, 255, 0.45);
    }}
    section {{
      padding: 10px 14px 8px 14px;
      margin: 8px 0 14px 0;
      background: #ffffff;
      border: 1px solid #dfe1df;
      border-radius: 10px;
      scroll-margin-top: 90px;
    }}
    .tabs-wrap {{
      position: sticky;
      top: 0;
      z-index: 1000;
      background: #d11947;
      border-radius: 0 0 8px 8px;
      padding: 8px 0 6px 0;
      margin-bottom: 14px;
      box-shadow: 0 2px 4px rgba(141, 0, 52, 0.18);
    }}
    .tabs {{
      display: flex;
      gap: 8px;
      align-items: center;
      flex-wrap: wrap;
      margin: 0 6px;
    }}
    .tabs a {{
      text-decoration: none;
      font-size: 13px;
      font-weight: 700;
      color: #ffffff;
      border: 1px solid rgba(255, 255, 255, 0.6);
      border-radius: 8px;
      padding: 6px 10px;
      background: rgba(255, 255, 255, 0.14);
    }}
    .tabs a:hover {{
      background: rgba(255, 255, 255, 0.26);
      border-color: #ffffff;
    }}
  </style>
</head>
<body>
<div class="top-banner">
<h1>SHARE-seq QC Report: {html.escape(sample)}</h1>
<p>Project path: <code>{html.escape(proj)}</code></p>
</div>
''')
    sample_parts.append('<div class="tabs-wrap"><div class="tabs">')
    for sec_id, sec_label in sample_tabs:
        sample_parts.append(f'<a href="#{html.escape(sec_id)}">{html.escape(sec_label)}</a>')
    sample_parts.append("</div></div>")
    sample_parts.append('<section id="sec-demux">')
    sample_parts.append("<h2>Demultiplexing</h2>")
    sample_parts.append(demux_stats_html_for_sample(demux_stats, sample))
    sample_parts.append(sample_links_block("FastQC HTML (demultiplexed)", sample_fastqc_demux))
    sample_parts.append(sample_links_block("FastQC HTML (trimmed)", sample_fastqc_trimmed))
    sample_parts.append("</section>")
    if has_atac:
        sample_parts.append('<section id="sec-atac">')
        sample_parts.append("<h2>ATAC QC</h2>")
        sample_parts.append(atac_pre_post_table(sample))
        sample_parts.append(sample_links_block("ATAC MultiQC report", atac_multiqc_report))
        sample_parts.append("</section>")
    if has_rna:
        sample_parts.append('<section id="sec-rna">')
        sample_parts.append("<h2>RNA QC</h2>")
        sample_parts.append(text_files_block("Log.final.out", sample_logs, max_lines=120))
        sample_parts.append(sample_image_block("Knee plots", sample_knee))
        sample_parts.append(text_files_block("Barcodes.stats", sample_barcodes, max_lines=120))
        sample_parts.append(table_files_block("GeneFull Summary.csv", sample_summary, max_rows=20))
        sample_parts.append(sample_image_block("Barnyard plots (if hybrid)", sample_barnyard))
        sample_parts.append("</section>")
    sample_parts.append("</body></html>")

    with open(os.path.join(sample_root, "index.html"), "w") as sf:
        sf.write("\n".join(sample_parts))

with zipfile.ZipFile("QC_Report_bundle.zip", "w", compression=zipfile.ZIP_DEFLATED) as zf:
    zf.write(out_path, arcname=out_path)
    for root, _, files in os.walk(assets_dir):
        for f in files:
            p = os.path.join(root, f)
            zf.write(p, arcname=p)
    for root, _, files in os.walk(per_sample_dir):
        for f in files:
            p = os.path.join(root, f)
            zf.write(p, arcname=p)
