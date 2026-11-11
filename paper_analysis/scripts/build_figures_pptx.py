"""Assemble the main manuscript figures (Fig. 1-2) into a single PowerPoint deck,
one A4-landscape slide per slide image, with the slide image plus the
corresponding slice of its legend text as a caption below.

Fig. 1 (5 panels, a-e) and Fig. 2 (6 panels, a-f) are each merges of what used
to be several separate figures. Per Simone, a figure's panels should sit
together on one composited image rather than just be labelled as parts of the
same figure across separate slides/pages -- so Fig. 1a-d and all of Fig. 2a-f
are each a single composite image (see scripts/compose_figures.py). Fig. 1e
(the cross-trait heatmap) is the one exception: its tall, many-row layout
doesn't composite cleanly alongside a-d, so it stays on its own slide.

Legend text here is a plain-text transcription of the legends in
MANUSCRIPT.tex's "Figure legends" section (kept in sync by hand -- there are
only 2 of them).

Run from the repo root (after scripts/compose_figures.py has produced the
composite images):
    python scripts/build_figures_pptx.py
"""

import os

from pptx import Presentation
from pptx.util import Cm, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor
from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG_DIR = os.path.join(ROOT, "figures")
OUT_PATH = os.path.join(FIG_DIR, "CREST-GV_main_figures.pptx")

PAGE_W, PAGE_H = Cm(29.7), Cm(21.0)  # A4 landscape
MARGIN = Cm(1.2)
TITLE_H = Cm(1.0)
LEGEND_H = Cm(4.2)
INK = RGBColor(0x1A, 0x1A, 0x1A)

FIGURES = [
    {
        "num": 1,
        "title": "Figure 1 (a-d)",
        "image": os.path.join(FIG_DIR, "Fig1_abcd_composite.png"),
        "legend": (
            "Fig. 1 | Overview and standardisation of CREST-GV enrichment scores "
            "(panel b is a placeholder figure). "
            "a, Enrichment framework. Trait-associated variants (GWAS lead variants and LD "
            "proxies, or a user-provided variant set) and background variants (1000 Genomes "
            "Project) are intersected with the open chromatin regions of each cell type in the "
            "reference atlas or a user-provided collection. For each cell type, the observed "
            "overlap proportion is compared with the background proportion to compute an "
            "enrichment Z-score, standardised within each data collection and converted to a "
            "right-tailed P-value. "
            "b, The CREST-GV web application, providing access to precomputed results for "
            "1,125 GWAS traits and analysis of user-provided variant sets and chromatin "
            "accessibility collections. "
            "c, Distributions of raw enrichment scores (Z) for each chromatin accessibility "
            "collection, showing collection-specific differences in scale, for red blood cell "
            "traits (one of the four traits highlighted in Results; chosen for consistency with "
            "panel e and Fig. 2, which use the same four traits). "
            "d, Distributions of standardised enrichment scores (Z̃) after within-collection "
            "standardisation, same trait and cell types. (Continues as Fig. 1e on the next slide.)"
        ),
        "placeholder": "Image not yet available -- run scripts/compose_figures.py first.",
    },
    {
        "num": 1,
        "title": "Figure 1 (e)",
        "image": os.path.join(FIG_DIR, "Fig3_cross_trait_heatmap.png"),
        "legend": (
            "e, CREST-GV enrichment across four GWAS traits: standardised CREST-GV "
            "enrichment scores (Z̃) for the top 10 cell types by score for each of heart rate, "
            "red blood cell traits, type 1 diabetes and type 2 diabetes (union across the four "
            "traits, not restricted to significance). Rows are grouped by the trait each cell "
            "type scores highest for; values in each cell give the standardised score of that "
            "cell type for every trait; asterisks mark Benjamini-Hochberg FDR q < 0.10 "
            "significance."
        ),
    },
    {
        "num": 2,
        "title": "Figure 2 (a-f)",
        "image": os.path.join(FIG_DIR, "Fig2_abcdef_composite.png"),
        "legend": (
            "Fig. 2 | Cell type rankings and analysis of user-provided data with CREST-GV. "
            "Top and bottom 10 standardised enrichment scores (Z̃) across all collections for "
            "a, heart rate; b, red blood cell traits; c, type 1 diabetes; and d, type 2 diabetes. "
            "Bars are coloured by sign; asterisks mark cell types with Benjamini-Hochberg "
            "FDR q < 0.10. "
            "e, Held-out collection re-upload: the cad collection re-analysed through the "
            "in-house-collection code path (same bigWig tracks and peak calls, processed as if "
            "newly uploaded), compared against its precomputed reference-atlas enrichment score "
            "for platelet_measurement. Dashed line, y = x; each point is one of the ten cad cell "
            "types. "
            "f, User-provided variant set: a GWAS trait not among the four traits highlighted in "
            "Results (osteoporosis), run fresh against all 11 reference collections; top 10 "
            "standardised scores shown."
        ),
        "placeholder": "Image not yet available -- run scripts/compose_figures.py first.",
    },
]


def add_title(slide, text):
    box = slide.shapes.add_textbox(MARGIN, Cm(0.3), PAGE_W - 2 * MARGIN, TITLE_H)
    tf = box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(22)
    p.font.bold = True
    p.font.color.rgb = INK
    return box


def add_image_fit(slide, image_path, area_left, area_top, area_w, area_h):
    with Image.open(image_path) as im:
        iw, ih = im.size
    area_ratio = area_w / area_h
    img_ratio = iw / ih
    if img_ratio > area_ratio:
        w = area_w
        h = int(area_w / img_ratio)
    else:
        h = area_h
        w = int(area_h * img_ratio)
    left = area_left + (area_w - w) // 2
    top = area_top + (area_h - h) // 2
    slide.shapes.add_picture(image_path, left, top, width=w, height=h)


def add_legend(slide, text):
    box = slide.shapes.add_textbox(MARGIN, PAGE_H - LEGEND_H - MARGIN, PAGE_W - 2 * MARGIN, LEGEND_H)
    tf = box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(11)
    p.font.color.rgb = INK
    p.alignment = PP_ALIGN.LEFT
    return box


def main():
    prs = Presentation()
    prs.slide_width = PAGE_W
    prs.slide_height = PAGE_H
    blank_layout = prs.slide_layouts[6]

    image_area_top = MARGIN + TITLE_H
    image_area_h = PAGE_H - TITLE_H - LEGEND_H - 2 * MARGIN
    image_area_w = PAGE_W - 2 * MARGIN

    for fig in FIGURES:
        slide = prs.slides.add_slide(blank_layout)
        add_title(slide, fig["title"])

        if fig["image"] and os.path.exists(fig["image"]):
            add_image_fit(slide, fig["image"], MARGIN, image_area_top, image_area_w, image_area_h)
        else:
            box = slide.shapes.add_textbox(MARGIN, image_area_top, image_area_w, image_area_h)
            tf = box.text_frame
            tf.word_wrap = True
            p = tf.paragraphs[0]
            p.text = fig.get("placeholder", "Image not available.")
            p.font.size = Pt(16)
            p.font.italic = True
            p.font.color.rgb = RGBColor(0x80, 0x80, 0x80)
            p.alignment = PP_ALIGN.CENTER

        add_legend(slide, fig["legend"])

    prs.save(OUT_PATH)
    print(f"wrote {OUT_PATH} ({len(FIGURES)} slides, A4 landscape)")


if __name__ == "__main__":
    main()
