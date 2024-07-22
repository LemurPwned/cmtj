import logging
import os
import re
import shutil

log = logging.getLogger("mkdocs")

KNOWN_REFLINK_MAP = {}


def on_startup(command, dirty, **kwargs):
    """Copy the README to the docs folder and replace reflinks."""
    shutil.copy("README.md", "./docs/index.md")
    shutil.copy("CHANGELOG.md", "./docs/changelog.md")
    if os.path.exists("./docs/assets"):
        shutil.rmtree("./docs/assets")
    shutil.copytree("assets", "./docs/assets")
    # replace all the reflinks
    with open("docs/index.md") as f:
        content = reflink_and_div_replace(f)
    with open("docs/index.md", "w") as f:
        f.write(content)


def reflink_and_div_replace(f):
    result = f.read()
    for reflink, new_reflink in KNOWN_REFLINK_MAP.items():
        print(f"Replacing {reflink} with {new_reflink}")
        result = result.replace(reflink, new_reflink)

    # find all the div tags
    divs = re.findall(r"<div.*?>", result)
    for div in divs:
        result = result.replace(div, "")
    # remove all the closing div tags
    divs = re.findall(r"</div>", result)
    for div in divs:
        result = result.replace(div, "")

    return result
