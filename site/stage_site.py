#!/usr/bin/env python3
"""Stage the existing cdid guide as a disposable pkgdown source tree.

Called by build_site.R; Python standard library only. Reads the source
checkout (never writes to it) and creates a new stage directory:

  pkg/DESCRIPTION, pkg/NAMESPACE, pkg/man/*.Rd   current source, unchanged
  pkg/NEWS.md                                   release news for pkgdown
  pkg/_pkgdown.yml                              copied from --config
  pkg/index.md                                  home page from docs/README.md
  pkg/vignettes/**/*.Rmd                        article adapters
  pkg/vignettes/figures/*.png                   existing plots, exact bytes
  site-manifest.json                            mapping and expected output

An adapter is the existing Markdown page with standard YAML front matter.
Only these mechanical changes are made: the page's H1 moves into the title,
relative links are mapped to the built HTML pages, a blank line is inserted
where pandoc needs one before a list, heading, table, quote or code fence,
and a heading gets an explicit {#id} where pandoc's automatic identifier
would differ from GitHub's. Code fences stay static code blocks; nothing is
evaluated. Any page, link, anchor or figure that cannot be mapped stops the
staging with an error before anything is written.

Usage: python3 -I stage_site.py --source DIR --stage NEW_DIR --config FILE --logo PNG
"""

import argparse
import hashlib
import html
import json
import os
import posixpath
import re
import shutil
import sys
import struct
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.parse import unquote, urlsplit

HOME = "docs/README.md"

# Source page (checkout relative) -> article name (articles/<name>.html).
ARTICLES = {
    "docs/get-started.md": "get-started",
    "docs/from-did.md": "from-did",
    "docs/examples/incomplete-panels.md": "examples/incomplete-panels",
    "docs/examples/identification.md": "examples/identification",
    "docs/examples/precision.md": "examples/precision",
    "docs/help/prepare-data.md": "help/prepare-data",
    "docs/help/investigate-missingness.md": "help/investigate-missingness",
    "docs/methods/treatment.md": "methods/treatment",
    "docs/examples/treatment-choice.md": "examples/treatment-choice",
    "docs/methods/missing-outcomes.md": "methods/missing-outcomes",
    "docs/methods/sampling-information.md": "methods/sampling-information",
    "docs/examples/observation-models.md": "examples/observation-models",
    "docs/examples/change-imputation.md": "examples/change-imputation",
    "docs/methods/parallel-trends.md": "methods/parallel-trends",
    "docs/examples/pt-all-comparison.md": "examples/pt-all-comparison",
    "docs/methods/gmm.md": "methods/gmm",
    "docs/examples/robustness.md": "examples/robustness",
    "docs/examples/inference.md": "examples/inference",
    "docs/reference/options.md": "reference/options",
    "docs/reference/results.md": "reference/results",
    "docs/help/troubleshooting.md": "help/troubleshooting",
    "validation/SETUP.md": "installation",
}

# Existing plots -> articles/figures/<file name>.
FIGURES = [
    "docs/outputs/clustered-cells.png",
    "docs/outputs/clustered-dynamic.png",
    "docs/outputs/get-started-cells.png",
    "docs/outputs/get-started-dynamic.png",
    "docs/outputs/incomplete-dr-dynamic.png",
]

# Linked files that are not pages: link to them in the private repository.
# Branch as in the Contributing link of docs/README.md.
REPO_BLOB = "https://github.com/dvbn/cdid/blob/main/"
REPO_FILES = [
    "notebooks/README.md",
    "notebooks/cdid_practical_guide.qmd",
    "notebooks/cdid_practical_guide.html",
    "notebooks/simulate_panel.R",
]

FENCE = re.compile(r"^\s*(`{3,}|~{3,})(.*)$")
CODE_SPAN = re.compile(r"(`+)(.+?)(?<!`)\1(?!`)")
LINK = re.compile(r"\]\((<?)([^)\s>]*)(>?)((?:\s+\"[^\"]*\")?)\)")
HEADING = re.compile(r"^(#{1,6})\s+(.*?)(?:\s+#+)?\s*$")
# Block starts that interrupt a paragraph on GitHub but not in pandoc.
LIST_START = re.compile(r"^([-*+]|1[.)])\s")
ANY_LIST = re.compile(r"^\s*([-*+]|\d+[.)])\s")


class StageError(Exception):
    pass


def site_path(rel):
    """Site location of a mapped source file, relative to the site root."""
    if rel in (HOME, "README.md"):
        return "index.html"
    if rel == "NEWS.md":
        return "news/index.html"
    if rel in ("man/figures/logo.png", "man/figures/logo.svg"):
        return posixpath.basename(rel)
    if rel in ARTICLES:
        return "articles/" + ARTICLES[rel] + ".html"
    if rel in FIGURES:
        return "articles/figures/" + posixpath.basename(rel)
    raise KeyError(rel)


def stage_path(rel):
    if rel == HOME:
        return "index.md"
    if rel == "NEWS.md":
        return "NEWS.md"
    if rel in ARTICLES:
        return "vignettes/" + ARTICLES[rel] + ".Rmd"
    return "vignettes/figures/" + posixpath.basename(rel)


def hashes(path):
    with open(path, "rb") as f:
        data = f.read()
    return {"bytes": len(data), "md5": hashlib.md5(data).hexdigest(),
            "sha256": hashlib.sha256(data).hexdigest()}


def is_external(target):
    return bool(re.match(r"[A-Za-z][A-Za-z0-9+.-]*:", target)) or target.startswith("//")


def mask_code(line):
    """Blank out inline code so that links inside it are left alone."""
    return CODE_SPAN.sub(lambda m: " " * len(m.group(0)), line)


def heading_text(content):
    s = re.sub(r"!?\[([^\]]*)\]\([^)]*\)", r"\1", content)
    s = re.sub(r"<[^>]+>", "", s)
    s = s.replace("`", "")
    s = re.sub(r"\*+", "", s)
    return s.strip()


def github_slug(text):
    s = "".join(c for c in text.lower() if c.isalnum() or c in " -_")
    return s.replace(" ", "-")


def pandoc_slug(text):
    s = "".join(c for c in text.lower() if c.isalnum() or c in "_-." or c.isspace())
    s = "-".join(s.split())
    while s and not s[0].isalpha():
        s = s[1:]
    return s or "section"


def dedupe(slug, seen):
    n = seen.get(slug, 0)
    seen[slug] = n + 1
    return slug if n == 0 else "%s-%d" % (slug, n)


def classify(text):
    """Return [(line, kind)] with kind 'text', 'open', 'code' or 'close'."""
    out, fence = [], None
    for line in text.split("\n"):
        m = FENCE.match(line)
        if fence:
            closes = m and m.group(1)[0] == fence[0] and len(m.group(1)) >= len(fence) \
                and not m.group(2).strip()
            out.append((line, "close" if closes else "code"))
            if closes:
                fence = None
        elif m:
            fence = m.group(1)
            out.append((line, "open"))
        else:
            out.append((line, "text"))
    if fence:
        raise StageError("unclosed code fence")
    return out


def split_title(rel, text):
    """Return (title, body). The page H1 becomes the title, once."""
    if rel in (HOME, "README.md", "NEWS.md"):
        return None, text
    first, _, rest = text.partition("\n")
    m = HEADING.match(first)
    if not m or m.group(1) != "#":
        raise StageError("first line is not the page H1")
    if rest.startswith("\n"):
        rest = rest[1:]
    return m.group(2), rest


def github_anchors(title, body):
    seen, anchors = {}, set()
    if title is not None:
        anchors.add(dedupe(github_slug(heading_text(title)), seen))
    for line, kind in classify(body):
        m = HEADING.match(line) if kind == "text" else None
        if m:
            anchors.add(dedupe(github_slug(heading_text(m.group(2))), seen))
    return anchors


def needs_blank(prev):
    """True if prev is a top-level paragraph line that pandoc would continue."""
    if not prev.strip() or prev[0] in " \t":
        return False
    return not (ANY_LIST.match(prev) or prev.startswith(("#", "|", ">", "<")))


def map_link(target, rel, anchors, record):
    if not target:
        raise StageError("empty link target")
    if is_external(target):
        return target
    url = urlsplit(html.unescape(target))
    path, frag = unquote(url.path), unquote(url.fragment)
    hash_ = "#" if url.fragment else ""
    if path.startswith("/") or "\\" in path:
        raise StageError("use a relative website link: " + target)
    if not path:
        dest_rel, new = rel, ""
    else:
        dest_rel = posixpath.normpath(posixpath.join(posixpath.dirname(rel), path))
        if dest_rel in REPO_FILES:
            if hash_:
                raise StageError("anchor on a repository file link: " + target)
            record["repo_links"].append({"target": dest_rel, "url": REPO_BLOB + dest_rel})
            return REPO_BLOB + dest_rel
        try:
            dest = site_path(dest_rel)
        except KeyError:
            raise StageError("unmapped relative link: " + target)
        if dest_rel in FIGURES:
            record["figures"].add(dest_rel)
        here = posixpath.dirname(site_path(rel)) or "."
        new = posixpath.relpath(dest, here)
        record["links_mapped"] += 1
    if hash_:
        if frag not in anchors.get(dest_rel, ()):
            raise StageError("anchor not found in %s: #%s" % (dest_rel, frag))
        record["anchors"].add((site_path(dest_rel), frag))
        new += "#" + frag
    if url.query:
        new = new.partition("#")[0] + "?" + url.query + ("#" + frag if hash_ else "")
    return new


def adapt(rel, title, body, anchors, record):
    """Return the staged Markdown body; record counts what was changed."""
    lines = classify(body)
    out, prev = [], ""
    gh_seen, pd_seen = {}, {}
    if title is not None:
        dedupe(github_slug(heading_text(title)), gh_seen)
    for i, (line, kind) in enumerate(lines):
        if kind in ("code", "close"):
            out.append(line)
            prev = line
            continue
        if kind == "open" and "{" in FENCE.match(line).group(2):
            raise StageError("knitr chunk header: " + line.strip())
        masked = mask_code(line)
        if re.search(r"(?<!`)`r[ #]", line):
            raise StageError("knitr inline code: " + line.strip())
        # Rebase raw HTML (including the logo in the home H1) just like Markdown.
        pieces, pos = [], 0
        for m in re.finditer(r"\b(?:href|src)\s*=\s*([\"'])(.*?)\1", masked, re.I):
            pieces.append(line[pos:m.start(2)])
            pieces.append(html.escape(map_link(m.group(2), rel, anchors, record), quote=True))
            pos = m.end(2)
        line = "".join(pieces) + line[pos:]
        masked = mask_code(line)
        ref = re.match(r"\s{0,3}\[[^\]]+\]:\s*(\S+)", masked)
        if ref and not is_external(ref.group(1)):
            raise StageError("relative reference link definition: " + ref.group(1))
        if kind == "text" and re.match(r"^(=+|-+)\s*$", line):
            if needs_blank(prev):
                raise StageError("setext heading after: " + prev.strip())
            following = lines[i + 1][0] if i + 1 < len(lines) else ""
            if line.strip() == "---" and following.strip():
                raise StageError("rule followed by text (pandoc reads YAML): " + following)
        heading = HEADING.match(line) if kind == "text" else None
        starts_block = kind == "open" or heading or LIST_START.match(line) \
            or line.startswith(("|", ">"))
        if starts_block and needs_blank(prev):
            out.append("")
            record["blank_lines_inserted"] += 1
        pieces, pos = [], 0
        for m in LINK.finditer(masked):
            pieces.append(line[pos:m.start(2)])
            pieces.append(map_link(m.group(2), rel, anchors, record))
            pos = m.end(2)
        line = "".join(pieces) + line[pos:]
        if heading and rel != "NEWS.md":
            text = heading_text(heading.group(2))
            gh = dedupe(github_slug(text), gh_seen)
            if gh != dedupe(pandoc_slug(text), pd_seen):
                m = HEADING.match(line)
                line = "%s %s {#%s}" % (m.group(1), m.group(2), gh)
                record["explicit_ids"] += 1
        out.append(line)
        prev = line
    return "\n".join(out)


def yaml_front_matter(title):
    # A JSON string is a valid YAML double-quoted scalar.
    return "---\ntitle: %s\n---\n\n" % json.dumps(title, ensure_ascii=False)


def check_complete(src):
    """Every Markdown page and plot under docs/ must be mapped."""
    found_md, found_png = set(), set()
    for root, dirs, files in os.walk(os.path.join(src, "docs")):
        for name in dirs:
            if not Path(root, name).resolve().is_relative_to(Path(src)):
                raise StageError("docs/ directory resolves outside --source")
        for name in files:
            rel = posixpath.relpath(os.path.join(root, name).replace(os.sep, "/"),
                                    src.replace(os.sep, "/"))
            if name.endswith(".md"):
                found_md.add(rel)
            elif name.endswith(".png"):
                found_png.add(rel)
    unmapped = sorted(found_md - set(ARTICLES) - {HOME}) + sorted(found_png - set(FIGURES))
    if unmapped:
        raise StageError("source files without a site mapping: " + ", ".join(unmapped))
    missing = [p for p in [HOME, "README.md", "NEWS.md"] + list(ARTICLES) + FIGURES
               if not os.path.isfile(os.path.join(src, p))]
    if missing:
        raise StageError("mapped source files missing: " + ", ".join(missing))


def check_config(config, rd_names):
    """The navbar and both indexes must list every article and help topic."""
    text = open(config, encoding="utf-8").read()
    if re.search(r"^(url|destination)\s*:", text, re.M):
        raise StageError("_pkgdown.yml must not set url or destination")
    problems = []
    for name in ARTICLES.values():
        if not re.search(r"^\s*-\s+%s\s*$" % re.escape(name), text, re.M):
            problems.append("articles index lacks " + name)
        if "articles/%s.html" % name not in text:
            problems.append("navbar lacks articles/%s.html" % name)
    for href in re.findall(r"articles/([\w/-]+)\.html", text):
        if href != "index" and href not in ARTICLES.values():
            problems.append("navbar links unknown article " + href)
    for name in rd_names:
        if not re.search(r"^\s*-\s+%s\s*$" % re.escape(name), text, re.M):
            problems.append("reference index lacks " + name)
    if problems:
        raise StageError("; ".join(problems))
    return text


PROVIDER = re.compile(r"\b(?:AI|OpenAI|Anthropic|Claude|Opus|ChatGPT|Gemini|Copilot|GPT(?:[- .]?\d[\w.-]*)?|agents?)\b", re.I)
INTERNAL_HEADING = re.compile(r"\b(?:CP\d{2,}\b|checkpoint\w*\b|development\s*:|development logs?\b|conventions for authors\b|decisions still open\b)", re.I)
INTERNAL_LABEL = re.compile(r"\bCP\d{2,}\b|\.runtime/(?:cp\d|rlib)", re.I)


def check_public(rel, text):
    if PROVIDER.search(text):
        raise StageError(rel + ": remove authoring-tool references before building the user website")
    if re.search(r"/(?:home|Users)/[^\s/]+/", text):
        raise StageError(rel + ": replace machine-specific paths before building the website")
    if re.search(r"\bfor authors of (?:the|this) guide\b", text, re.I):
        raise StageError(rel + ": replace authoring notes with user documentation before building")
    for line, kind in classify(text):
        if kind == "text" and HEADING.match(line) and INTERNAL_HEADING.search(line):
            raise StageError(rel + ": replace development headings with user documentation before building")
        if kind == "text" and INTERNAL_LABEL.search(line):
            raise StageError(rel + ": replace development labels with user documentation before building")


def source_file(src, rel):
    path = Path(src, rel)
    if not path.resolve().is_relative_to(Path(src)):
        raise StageError("source file resolves outside --source: " + rel)
    if not path.is_file():
        raise StageError("source file missing: " + rel)
    return path


def check_news(text, version):
    headings = [m for line, kind in classify(text)
                if kind == "text" and (m := HEADING.match(line))]
    if not headings or len(headings[0].group(1)) not in (1, 2):
        raise StageError("NEWS.md: start with a release heading")
    level = len(headings[0].group(1))
    if any(len(m.group(1)) < level for m in headings):
        raise StageError("NEWS.md: use consistent release heading levels")
    headings = [m.group(2) for m in headings if len(m.group(1)) == level]
    versions = []
    for heading in headings:
        m = re.fullmatch(r"(?:cdid\s+)?(\d+\.\d+\.\d+(?:\.\d+)?)(?:\s+\(\d{4}-\d{2}-\d{2}\))?", heading)
        if not m:
            raise StageError("NEWS.md: use release headings such as '# cdid " + version + "'")
        versions.append(m.group(1))
    if not versions or versions[0] != version or len(set(versions)) != len(versions):
        raise StageError("NEWS.md: the first release heading must match DESCRIPTION Version, with no duplicate versions")


def new_stage(src, stage):
    path = Path(stage).absolute()
    if path.exists() or path.is_symlink():
        raise StageError("stage directory already exists: " + str(path))
    parent = path.parent.resolve(strict=True)
    path = parent / path.name
    if path.is_relative_to(Path(src)):
        raise StageError("--stage must be outside --source")
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--source", required=True)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--logo", required=True, help="PNG supplied for the website")
    args = ap.parse_args()
    src = os.path.realpath(args.source)
    if not os.path.isdir(src):
        raise StageError("--source must be a directory")
    stage = str(new_stage(src, args.stage))
    for rel in ("docs", "man"):
        if not Path(src, rel).resolve().is_relative_to(Path(src)):
            raise StageError(rel + ": directory resolves outside --source")

    with source_file(src, "DESCRIPTION").open(encoding="utf-8") as f:
        description = f.read()
        if not re.search(r"^Package:\s*cdid\s*$", description, re.M):
            raise StageError("--source is not the cdid checkout")
    version = re.search(r"^Version:\s*(\S+)\s*$", description, re.M)
    if not version or not re.fullmatch(r"\d+\.\d+\.\d+(?:\.\d+)?", version.group(1)):
        raise StageError("DESCRIPTION: supply a release Version before building")
    version = version.group(1)
    if int(version.split(".")[-1]) >= 9000:
        raise StageError("DESCRIPTION: supply a release Version before building")
    check_public("DESCRIPTION", description)
    news = source_file(src, "NEWS.md").read_text(encoding="utf-8")
    check_news(news, version)
    check_public("NEWS.md", news)
    # Check symlinks before any mapped document is opened.
    for rel in [HOME, "README.md", "NAMESPACE"] + list(ARTICLES) + FIGURES:
        source_file(src, rel)
    check_complete(src)

    rd_files = sorted(n for n in os.listdir(os.path.join(src, "man")) if n.endswith(".Rd"))
    rd_names = []
    for name in rd_files:
        text = source_file(src, "man/" + name).read_text(encoding="utf-8")
        if "\\Sexpr" in text:
            raise StageError("man/%s uses \\Sexpr, which would run R code" % name)
        check_public("man/" + name, text)
        if re.search(r"\\(?:title|section)\{[^}]*" + INTERNAL_HEADING.pattern, text, re.I):
            raise StageError("man/" + name + ": replace development headings with user documentation before building")
        topic = re.search(r"\\name\{([^}]+)\}", text)
        if not topic:
            raise StageError("man/" + name + ": missing help topic name")
        rd_names.append(topic.group(1))
    check_config(args.config, rd_names)

    # Convert everything in memory first, so a failure writes nothing.
    pages = [HOME] + list(ARTICLES) + ["NEWS.md"]
    split, anchors = {}, {}
    for rel in pages:
        text = open(os.path.join(src, rel), encoding="utf-8").read()
        try:
            check_public(rel, text)
            split[rel] = split_title(rel, text)
            anchors[rel] = github_anchors(*split[rel])
        except StageError as e:
            raise StageError("%s: %s" % (rel, e))
    anchors["README.md"] = github_anchors(*split_title("README.md", source_file(src, "README.md").read_text(encoding="utf-8")))
    staged, page_records = {}, []
    all_anchors, figures_used, repo_links = set(), set(), []
    for rel in pages:
        title, body = split[rel]
        record = {"links_mapped": 0, "blank_lines_inserted": 0, "explicit_ids": 0,
                  "anchors": set(), "figures": set(), "repo_links": []}
        try:
            text = adapt(rel, title, body, anchors, record)
        except StageError as e:
            raise StageError("%s: %s" % (rel, e))
        title_anchor = '<a id="%s"></a>\n\n' % html.escape(github_slug(heading_text(title)), quote=True) if title else ""
        staged[rel] = text if rel in (HOME, "NEWS.md") else yaml_front_matter(title) + title_anchor + text
        all_anchors |= record["anchors"]
        figures_used |= record["figures"]
        repo_links += [dict(r, source=rel) for r in record["repo_links"]]
        page_records.append({
            "source": rel, "stage": "pkg/" + stage_path(rel), "site": site_path(rel),
            "title": title, "links_mapped": record["links_mapped"],
            "blank_lines_inserted": record["blank_lines_inserted"],
            "explicit_ids": record["explicit_ids"]})
    unused = sorted(set(FIGURES) - figures_used)
    if unused:
        raise StageError("figures not referenced by any page: " + ", ".join(unused))
    # The two home sources must agree after rebasing their relative URLs.
    root_body = source_file(src, "README.md").read_text(encoding="utf-8")
    check_public("README.md", root_body)
    root_record = {"links_mapped": 0, "blank_lines_inserted": 0, "explicit_ids": 0,
                   "anchors": set(), "figures": set(), "repo_links": []}
    if adapt("README.md", None, root_body, anchors, root_record) != staged[HOME]:
        raise StageError("README.md and docs/README.md must have the same content apart from relative link paths")
    missing_repo = sorted({r["target"] for r in repo_links
                           if not os.path.isfile(os.path.join(src, r["target"]))})

    for rel in FIGURES:
        with open(os.path.join(src, rel), "rb") as f:
            if f.read(8) != b"\x89PNG\r\n\x1a\n":
                raise StageError(rel + " is not a PNG file")
    logo = Path(args.logo).resolve(strict=True)
    logo_bytes = logo.read_bytes()
    if len(logo_bytes) < 24 or logo_bytes[:8] != b"\x89PNG\r\n\x1a\n" or logo_bytes[12:16] != b"IHDR":
        raise StageError("--logo is not a PNG file")
    if not all(struct.unpack(">II", logo_bytes[16:24])):
        raise StageError("--logo has invalid PNG dimensions")
    logos = [(logo, "logo.png")]
    svg = Path(args.logo).absolute().with_suffix(".svg")
    if svg.is_symlink():
        raise StageError("the companion SVG must be an ordinary file")
    if svg.exists():
        svg_text = svg.read_text(encoding="utf-8")
        if re.search(r"<!DOCTYPE|<!ENTITY", svg_text, re.I):
            raise StageError("the companion SVG must be self-contained")
        try:
            tree = ET.fromstring(svg_text)
        except ET.ParseError as error:
            raise StageError("invalid companion SVG: " + str(error))
        if tree.tag.split("}")[-1] != "svg":
            raise StageError("the companion SVG must have an svg root")
        for element in tree.iter():
            if element.tag.split("}")[-1] in ("script", "foreignObject"):
                raise StageError("the companion SVG must be a static image")
            for key, value in element.attrib.items():
                if key.lower().startswith("on"):
                    raise StageError("the companion SVG must be a static image")
                if key.split("}")[-1] in ("href", "src") and not value.startswith(("#", "data:")):
                    raise StageError("the companion SVG must not reference external files")
        for target in re.findall(r"url\(\s*[\"']?([^\"'\s)]+)", svg_text, re.I):
            if not target.startswith(("#", "data:")):
                raise StageError("the companion SVG must not reference external files")
        if re.search(r"@import\b", svg_text, re.I):
            raise StageError("the companion SVG must not import stylesheets")
        check_public("companion SVG", " ".join(tree.itertext()))
        logos.append((svg, "logo.svg"))

    # Write the stage. The parent directory must already exist.
    os.mkdir(stage)
    pkg = os.path.join(stage, "pkg")
    os.makedirs(os.path.join(pkg, "man"))
    for name in ["DESCRIPTION", "NAMESPACE"]:
        shutil.copyfile(os.path.join(src, name), os.path.join(pkg, name))
    for name in rd_files:
        shutil.copyfile(os.path.join(src, "man", name), os.path.join(pkg, "man", name))
    os.makedirs(os.path.join(pkg, "man", "figures"))
    for logo_file, name in logos:
        shutil.copyfile(logo_file, os.path.join(pkg, "man", "figures", name))
        # pkgdown copies root README assets separately from its template logo.
        shutil.copyfile(logo_file, os.path.join(pkg, name))
    shutil.copyfile(args.config, os.path.join(pkg, "_pkgdown.yml"))
    for rel, text in staged.items():
        dst = os.path.join(pkg, stage_path(rel))
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        with open(dst, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
    for rel in FIGURES:
        dst = os.path.join(pkg, stage_path(rel))
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(os.path.join(src, rel), dst)

    read = [HOME, "README.md", "NEWS.md"] + list(ARTICLES) + FIGURES + ["DESCRIPTION", "NAMESPACE"] + \
        ["man/" + n for n in rd_files]
    reference = ["reference/%s.html" % n[:-3] for n in rd_files]
    manifest = {
        "source_files": [dict(path=p, **hashes(os.path.join(src, p))) for p in read],
        "pages": page_records,
        "figures": [dict(source=rel, stage="pkg/" + stage_path(rel), site=site_path(rel),
                         **hashes(os.path.join(src, rel))) for rel in FIGURES],
        "reference_pages": reference,
        "reference_topics": rd_names,
        "version": version,
        "logo": hashes(logo),
        "logos": [dict(site=name, **hashes(path)) for path, name in logos],
        "anchors": [{"site": s, "id": a} for s, a in sorted(all_anchors)],
        "repo_links": repo_links,
        "repo_files_missing_from_source": missing_repo,
        "expected_site_files": ["index.html", "authors.html", "articles/index.html",
                                "reference/index.html", "news/index.html", "search.json"] +
                               ["articles/%s.html" % a for a in ARTICLES.values()] + reference,
    }
    with open(os.path.join(stage, "site-manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
        f.write("\n")
    for target in missing_repo:
        print("Warning: linked file not in the source checkout: " + target, file=sys.stderr)
    print("Staged %d pages, %d figures and %d help topics." %
          (len(pages), len(FIGURES), len(rd_files)))


if __name__ == "__main__":
    try:
        main()
    except (StageError, OSError) as e:
        sys.exit("Staging failed: %s" % e)
