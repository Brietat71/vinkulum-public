#!/usr/bin/env python3
"""Rend le rapport scientifique canonique en PDF, sans accès réseau.

Usage : python ci/rend_recherche.py
Dépendance de rendu facultative : reportlab==5.0.1.
Polices : DejaVu Sans, paquet système fonts-dejavu-core.
La source et le registre de provenance restent indépendants du moteur Vinkulum.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
from pathlib import Path
import re

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    HRFlowable,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/recherche/2026-09-07/report-source.md"
OUTPUT = ROOT / "docs/ETUDE_SCIENTIFIQUE_VINKULUM_2001_2026.pdf"
INK = colors.HexColor("#192D40")
ACCENT = colors.HexColor("#086E78")
LIGHT = colors.HexColor("#EFF5F7")
MUTED = colors.HexColor("#536576")
RULE = colors.HexColor("#D0DCE2")


def fonts() -> None:
    directory = Path("/usr/share/fonts/truetype/dejavu")
    for name, filename in (
        ("DejaVu", "DejaVuSans.ttf"),
        ("DejaVu-Bold", "DejaVuSans-Bold.ttf"),
        ("DejaVu-Italic", "DejaVuSans-Oblique.ttf"),
        ("DejaVu-BoldItalic", "DejaVuSans-BoldOblique.ttf"),
    ):
        pdfmetrics.registerFont(TTFont(name, str(directory / filename)))
    pdfmetrics.registerFontFamily(
        "DejaVu",
        normal="DejaVu",
        bold="DejaVu-Bold",
        italic="DejaVu-Italic",
        boldItalic="DejaVu-BoldItalic",
    )


def inline(value: str) -> str:
    """Le sous-ensemble Markdown du rapport : liens, gras et code simple."""
    result = []
    pattern = re.compile(r"\[([^\]]+)\]\((https?://[^\s]+)\)|\*\*(.+?)\*\*|`([^`]+)`")
    pos = 0
    for match in pattern.finditer(value):
        result.append(html.escape(value[pos : match.start()]))
        label, url, bold, code = match.groups()
        if url:
            result.append(
                '<a href="{}" color="#086E78"><u>{}</u></a>'.format(
                    html.escape(url, quote=True), html.escape(label)
                )
            )
        elif bold:
            result.append(f"<b>{html.escape(bold)}</b>")
        else:
            result.append(f"<font color='#086E78'>{html.escape(code)}</font>")
        pos = match.end()
    result.append(html.escape(value[pos:]))
    return "".join(result)


def styles() -> dict[str, ParagraphStyle]:
    body = ParagraphStyle(
        "Body",
        fontName="DejaVu",
        fontSize=10,
        leading=14.4,
        textColor=INK,
        spaceAfter=8,
        allowWidows=0,
        allowOrphans=0,
        alignment=TA_LEFT,
    )
    return {
        "body": body,
        "title": ParagraphStyle(
            "Title", parent=body, fontName="DejaVu-Bold", fontSize=22,
            leading=27, spaceAfter=13, keepWithNext=True,
        ),
        "h2": ParagraphStyle(
            "Section", parent=body, fontName="DejaVu-Bold", fontSize=17,
            leading=22, spaceBefore=4, spaceAfter=14, keepWithNext=True,
        ),
        "subtitle": ParagraphStyle(
            "Subtitle", parent=body, fontSize=12, leading=17,
            textColor=ACCENT, spaceAfter=14, keepWithNext=True,
        ),
        "h3": ParagraphStyle(
            "Subsection", parent=body, fontName="DejaVu-Bold", fontSize=11,
            leading=15, textColor=ACCENT, spaceBefore=10, spaceAfter=6,
            keepWithNext=True,
        ),
        "cell": ParagraphStyle(
            "Cell", parent=body, fontSize=8.8, leading=12.2, spaceAfter=0,
        ),
        "thead": ParagraphStyle(
            "TableHead", parent=body, fontName="DejaVu-Bold", fontSize=8.7,
            leading=12, spaceAfter=0, textColor=colors.white,
        ),
        "list": ParagraphStyle(
            "List", parent=body, leftIndent=15, firstLineIndent=-15,
            spaceAfter=7,
        ),
        "reference": ParagraphStyle(
            "Reference", parent=body, fontSize=8.4, leading=11.6,
            spaceAfter=9, allowWidows=0, allowOrphans=0,
        ),
    }


class Report(BaseDocTemplate):
    def __init__(self, filename: Path, title: str | None = None):
        super().__init__(
            str(filename),
            pagesize=A4,
            leftMargin=47,
            rightMargin=47,
            topMargin=51,
            bottomMargin=45,
            pageCompression=1,
            title=title or "Vinkulum : avancées scientifiques 2001–2026",
            author="Projet Vinkulum",
            subject="Recherche scientifique et décisions d’architecture — 7 septembre 2026",
            creator="Vinkulum — ci/rend_recherche.py",
            invariant=1,
        )
        frame = Frame(
            self.leftMargin, self.bottomMargin, self.width, self.height,
            leftPadding=0, bottomPadding=0, rightPadding=0, topPadding=0,
        )
        self.addPageTemplates(PageTemplate(id="Report", frames=frame, onPage=self.decorate))
        self.heading_records = []

    def decorate(self, canvas, doc) -> None:
        # Date éditoriale fixe : le rendu reproductible ne doit pas afficher 2000.
        canvas.setDateFormatter(lambda *args: "D:20260907000000+00'00'")
        canvas.saveState()
        width, height = A4
        canvas.setFillColor(ACCENT)
        canvas.rect(0, height - 6, width, 6, fill=1, stroke=0)
        canvas.setFont("DejaVu-Bold", 8)
        canvas.setFillColor(INK)
        canvas.drawString(self.leftMargin, height - 29, "VINKULUM")
        canvas.setFont("DejaVu", 8)
        canvas.setFillColor(MUTED)
        canvas.drawRightString(width - self.rightMargin, height - 29, "RECHERCHE SCIENTIFIQUE · 2001–2026")
        canvas.setStrokeColor(RULE)
        canvas.setLineWidth(0.5)
        canvas.line(self.leftMargin, 33, width - self.rightMargin, 33)
        canvas.setFont("DejaVu", 7.5)
        canvas.drawString(self.leftMargin, 21, "État des connaissances et du projet au 7 septembre 2026")
        canvas.drawRightString(width - self.rightMargin, 21, str(doc.page))
        canvas.restoreState()

    def afterFlowable(self, flowable) -> None:
        if not isinstance(flowable, Paragraph):
            return
        if flowable.style.name not in {"Title", "Section"}:
            return
        label = flowable.getPlainText()
        key = f"section-{len(self.heading_records)}"
        self.canv.bookmarkPage(key)
        self.canv.addOutlineEntry(label, key, level=0, closed=False)
        self.heading_records.append({"title": label, "page": self.page})


def table(lines: list[str], width: float, style: dict) -> Table:
    rows = [[s.strip() for s in line.strip().strip("|").split("|")] for line in lines]
    rows = [row for row in rows if not all(re.fullmatch(r":?-+:?", s) for s in row)]
    columns = len(rows[0])
    if any(len(row) != columns for row in rows):
        raise ValueError("Tableau Markdown non rectangulaire")
    ratios = [0.27, 0.39, 0.34] if columns == 3 else [1 / columns] * columns
    if columns == 3 and rows[0][0] == "Références":
        ratios = [0.21, 0.47, 0.32]
    data = [
        [Paragraph(inline(cell), style["thead" if i == 0 else "cell"]) for cell in row]
        for i, row in enumerate(rows)
    ]
    item = Table(data, colWidths=[r * width for r in ratios], repeatRows=1, hAlign="LEFT")
    item.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), INK),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [LIGHT, colors.white]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 7),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
        ("LINEBELOW", (0, -1), (-1, -1), 0.5, RULE),
    ]))
    return item


def build(source: Path, output: Path, title: str | None = None) -> dict:
    fonts()
    style = styles()
    doc = Report(output, title)
    lines = source.read_text(encoding="utf-8").splitlines()
    story = []
    i = 0
    references = False
    content = "\n".join(lines)
    if "{{" in content or re.search(r"\bturn\d+(?:view|search|fetch)\d+\b", content):
        raise ValueError("Marqueur interne ou substitution inachevée dans la source")
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line:
            continue
        if line == "---":
            story.append(PageBreak())
            continue
        if line.startswith("|"):
            block = [line]
            while i < len(lines) and lines[i].strip().startswith("|"):
                block.append(lines[i].strip())
                i += 1
            story.extend([Spacer(1, 3), table(block, doc.width, style), Spacer(1, 11)])
            continue
        if line.startswith("# "):
            story.append(Paragraph(inline(line[2:]), style["title"]))
            continue
        if line.startswith("## "):
            label = line[3:]
            references = label == "Sources et accès"
            kind = "subtitle" if label.startswith("Recherche et décisions") else "h2"
            story.append(Paragraph(inline(label), style[kind]))
            if kind == "h2":
                story.append(HRFlowable(width="100%", thickness=0.5, color=RULE, spaceAfter=12))
            continue
        if line.startswith("### "):
            story.append(Paragraph(inline(line[4:]), style["h3"]))
            continue
        block = [line]
        while i < len(lines) and lines[i].strip():
            if re.match(r"^(#|\||---$|\d+\. )", lines[i].strip()):
                break
            block.append(lines[i].strip())
            i += 1
        value = " ".join(block)
        kind = "reference" if references else "list" if re.match(r"^\d+\. ", value) else "body"
        story.append(Paragraph(inline(value), style[kind]))
    output.parent.mkdir(parents=True, exist_ok=True)
    doc.build(story)
    return {
        "source": str(source.relative_to(ROOT)) if source.is_relative_to(ROOT) else str(source),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "pdf": str(output.relative_to(ROOT)) if output.is_relative_to(ROOT) else str(output),
        "pdf_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "bytes": output.stat().st_size,
        "headings": doc.heading_records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--sortie", type=Path, default=OUTPUT)
    parser.add_argument("--titre", help="Titre des métadonnées PDF")
    args = parser.parse_args()
    print(json.dumps(build(args.source.resolve(), args.sortie.resolve(), args.titre), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
