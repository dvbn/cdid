#!/usr/bin/env python3
"""Check a rendered local pkgdown site without starting a browser or server."""

import argparse
import hashlib
import json
import re
import runpy
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

# -I intentionally omits the script directory from Python's import path.
_stage = runpy.run_path(str(Path(__file__).with_name("stage_site.py")))
INTERNAL_HEADING, INTERNAL_LABEL, PROVIDER, StageError = (_stage[k] for k in ("INTERNAL_HEADING", "INTERNAL_LABEL", "PROVIDER", "StageError"))

CSS_URL = re.compile(r"url\(\s*['\"]?([^'\"\s)]+)['\"]?\s*\)|@import\s+['\"]([^'\"]+)['\"]", re.I)


class Page(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.ids, self.refs, self.text, self.prose, self.headings = set(), [], [], [], []
        self.heading, self.hidden = None, 0
        self.pre_depth = 0
        self.noindex = False

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.ids.add(attrs["id"])
        if tag == "a" and attrs.get("name"):
            self.ids.add(attrs["name"])
        if tag in ("script", "style"):
            self.hidden += 1
        if tag == "pre":
            self.pre_depth += 1
        if re.fullmatch(r"h[1-6]", tag):
            self.heading = []
        if tag == "meta" and attrs.get("name", "").lower() == "robots":
            self.noindex = "noindex" in attrs.get("content", "").lower()
        for attr in ("href", "src", "poster", "data-search-index", "xlink:href"):
            if attrs.get(attr):
                runtime = tag != "a" or attr != "href"
                self.refs.append((attrs[attr], runtime))
        for target in re.findall(r"(?:^|,)\s*(\S+)(?:\s+[\d.]+[wx])?", attrs.get("srcset", "")):
            self.refs.append((target, True))
        if attrs.get("data-bs-target", "").startswith("#"):
            self.refs.append((attrs["data-bs-target"], False))
        self.refs += [(a or b, True) for a, b in CSS_URL.findall(attrs.get("style", ""))]

    def handle_endtag(self, tag):
        if tag in ("script", "style"):
            self.hidden = max(0, self.hidden - 1)
        if tag == "pre":
            self.pre_depth = max(0, self.pre_depth - 1)
        if re.fullmatch(r"h[1-6]", tag) and self.heading is not None:
            self.headings.append("".join(self.heading))
            self.heading = None

    def handle_data(self, text):
        if not self.hidden:
            self.text.append(text)
            if not self.pre_depth:
                self.prose.append(text)
        if self.heading is not None:
            self.heading.append(text)


def local_target(root, current, target, runtime=False):
    url = urlsplit(target)
    if url.scheme or url.netloc:
        if url.scheme.lower() == "data":
            return None
        if runtime or url.scheme.lower() not in ("https", "http", "mailto", "tel"):
            raise StageError(current + ": nonlocal asset or unsupported link: " + target)
        return None
    path = unquote(url.path)
    if "\\" in path:
        raise StageError(current + ": invalid local URL: " + target)
    dest = root / path.lstrip("/") if path.startswith("/") else root / Path(current).parent / path
    if not path:
        dest = root / current
    dest = dest.resolve()
    if not dest.is_relative_to(root):
        raise StageError(current + ": link escapes the website: " + target)
    if dest.is_dir():
        dest /= "index.html"
    if not dest.is_file():
        raise StageError(current + ": local target missing: " + target)
    if dest.suffix.lower() in (".md", ".rmd", ".rd"):
        raise StageError(current + ": link still points to source documentation: " + target)
    return dest.relative_to(root).as_posix(), unquote(url.fragment)


def check_site(output, manifest, host_paths=()):
    root = Path(output).resolve(strict=True)
    files = sorted(root.rglob("*"))
    if any(p.is_symlink() for p in files):
        raise StageError("website output must not contain symbolic links")
    for rel in manifest["expected_site_files"]:
        if not (root / rel).is_file():
            raise StageError("expected website file missing: " + rel)
    for entry in manifest["figures"] + manifest["logos"]:
        p = root / entry["site"]
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest() != entry["sha256"]:
            raise StageError("asset differs from its supplied source: " + entry["site"])

    pages, svg_ids, refs = {}, {}, []
    text_suffixes = {".html", ".css", ".js", ".json", ".txt", ".xml", ".svg", ".yml", ".map", ".webmanifest"}
    for path in files:
        if not path.is_file() or path.suffix.lower() not in text_suffixes:
            continue
        current = path.relative_to(root).as_posix()
        text = path.read_text(encoding="utf-8")
        # Absolute paths may also occur in JSON-escaped form or URL-encoded HTML.
        decoded = unquote(text.replace("\\/", "/"))
        for host in host_paths:
            if len(host) > 1 and host in decoded:
                raise StageError(current + ": contains a machine-specific source path")
        if re.search(r"/(?:home|Users)/[^\s/\"'<>]+/|\.site-build/|file://", decoded):
            raise StageError(current + ": contains a machine-specific or staging path")
        if path.suffix == ".html":
            page = Page()
            page.feed(text)
            page.close()
            pages[current] = page
            if not page.noindex:
                raise StageError(current + ": missing noindex metadata")
            if PROVIDER.search(" ".join(page.text)):
                raise StageError(current + ": contains authoring-tool references")
            if any(INTERNAL_HEADING.search(h) for h in page.headings):
                raise StageError(current + ": contains a development heading")
            if INTERNAL_LABEL.search(" ".join(page.prose)):
                raise StageError(current + ": contains development labels outside recorded output")
            if re.search(r"(?:github\.com/[^\s\"<>]+/blob/[^\s\"<>]+/(?:vignettes/|\.site-build/))", text):
                raise StageError(current + ": source link points to a generated adapter")
            refs += [(current, target, runtime) for target, runtime in page.refs]
        elif path.suffix == ".css":
            refs += [(current, a or b, True) for a, b in CSS_URL.findall(re.sub(r"/\*.*?\*/", "", text, flags=re.S))]
        elif path.suffix == ".webmanifest":
            data = json.loads(text)
            refs += [(current, icon["src"], True) for icon in data.get("icons", [])]
        elif path.suffix == ".svg":
            svg = Page()
            svg.feed(text)
            svg_ids[current] = svg.ids
            if PROVIDER.search(" ".join(svg.text)):
                raise StageError(current + ": contains authoring-tool references")
            refs += [(current, target, True) for target, _ in svg.refs]
            refs += [(current, a or b, True) for a, b in CSS_URL.findall(text)]

    for anchor in manifest["anchors"]:
        if anchor["id"] not in pages[anchor["site"]].ids:
            raise StageError(anchor["site"] + ": missing source anchor #" + anchor["id"])
    search = json.loads((root / "search.json").read_text(encoding="utf-8"))
    if not isinstance(search, list) or not search:
        raise StageError("search.json must contain a nonempty local search index")
    indexed = set()
    for row in search:
        # pkgdown 2.2.0 serializes empty sections as a path-only empty array.
        # They have no searchable fields; real entries still require local targets.
        if row == {"path": []}:
            continue
        if not isinstance(row, dict) or not isinstance(row.get("path"), str):
            raise StageError("search.json has an unexpected entry format")
        target = local_target(root, "index.html", row["path"], True)
        if target is None:
            raise StageError("search.json must use local page paths")
        indexed.add(target[0])
        fragment = row.get("id", "")
        refs.append(("index.html", row["path"] + ("#" + fragment if fragment else ""), True))
    expected_index = {p["site"] for p in manifest["pages"]} | set(manifest["reference_pages"])
    if not expected_index <= indexed:
        raise StageError("search.json omits pages: " + ", ".join(sorted(expected_index - indexed)))

    checked = 0
    for current, target, runtime in refs:
        resolved = local_target(root, current, target, runtime)
        if resolved is None:
            continue
        dest, fragment = resolved
        if fragment and dest in pages and fragment not in pages[dest].ids:
            raise StageError(current + ": anchor missing: " + target)
        if fragment and dest in svg_ids and fragment not in svg_ids[dest]:
            raise StageError(current + ": SVG anchor missing: " + target)
        checked += 1
    print("Checked %d HTML pages, %d local references, recorded plots, logo and search index. Browser review remains." %
          (len(pages), checked))
    return {"html_pages": len(pages), "local_references": checked, "search_pages": len(indexed)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--host-path", action="append", default=[])
    args = parser.parse_args()
    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    check_site(args.output, manifest, args.host_path)


if __name__ == "__main__":
    try:
        main()
    except (StageError, OSError, ValueError) as error:
        sys.exit("Website check failed: " + str(error))
