#!/usr/bin/env python3
"""Shared cell barcode normalization for SHARE-seq post-pipeline scripts."""

from __future__ import annotations

import re
from typing import Optional


def normalize_barcode(bc: str, length: int = 24) -> Optional[str]:
    if bc is None:
        return None
    s = str(bc).strip().upper()
    if not s:
        return None
    # ArchR / STAR may append sample suffix after '#'
    if "#" in s:
        s = s.split("#")[-1]
    s = re.sub(r"[^ACGTN]", "", s)
    if len(s) < length:
        return None
    if len(s) > length:
        s = s[:length]
    return s
