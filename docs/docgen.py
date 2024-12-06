import glob
import logging
import os
import re
from dataclasses import dataclass

log = logging.getLogger("mkdocs")

py_signature = r"(def (.+?) -> ([\'\[\]\,\sA-z]+)\:)"
c_py_signature = re.compile(py_signature)

# old regex that does not allow for optional docstrings and quoted rtype
# joint_py = r"(?s)(def (.+?) -> ([\"\'\[\]\,\sA-z]+)\:)\n{0,}\s{8}(\"{3}(.+?)\"{3})"
joint_py = r"(?s)(def (.+?) -> ([\"\'\[\]\,\sA-z]+)\:)(?:\n\s*(\"{3}(.+?)\"{3}))?"
c_joint_py = re.compile(joint_py)

pydoc_regex = r"(?s)(\"{3}(.+?)\"{3})"
c_pydoc_rgx = re.compile(pydoc_regex)
arg_pydoc = r"(\:param ([A-z0-9]+)\:(.+)\n)"
py_arg_rgx = re.compile(arg_pydoc)

cdoc_regex = r"(?s)(\/\*{2}(.+?)\*\/{1})"  # ?s is an inline DOTALL flag
arg_cdoc = r"(\@param ([A-z0-9]+)\:(.+)\n)"
c_cdoc_rgx = re.compile(cdoc_regex)
c_arg_rgx = re.compile(arg_cdoc)

GEN_FOLDER = "gen-docs"


@dataclass
class PythonDocstring:
    signature: str
    rtype: str
    docstring: str

    def extract_signature_types(self):
        type_map = {}
        rtype_map = {}
        first_bracket = self.signature.index("(")
        second_bracket = self.signature.index(")")
        args = self.signature[first_bracket + 1 : second_bracket].split(",")
        args = [arg.strip() for arg in args if arg != "self"]
        for arg in args:
            if ":" in arg:
                key_, type_ = arg.split(":")
                key_ = key_.strip()
                type_ = type_.strip()
                default_ = "-"
                if "=" in type_:
                    type_, default_ = type_.split("=")
                    type_ = type_.strip()
                    default_ = default_.strip()
                type_map[key_] = type_
                rtype_map[key_] = default_
        return type_map, rtype_map

    def py_signature_to_markdown(self):
        # form type map first
        type_map, rtype_map = self.extract_signature_types()

        arg_template = "**`{}`** | `{}` | {} | `{}`"
        table = (
            """Name | Type | Description | Default\n"""
            """------ | ---- | ----------- | -------"""
        )
        table = "#### **Parameters** \n" + table
        arg_count = 0
        for arg in py_arg_rgx.findall(self.docstring):
            if arg:
                arg_count += 1
                real_arg = arg[1].replace("\n", "")
                arg_desc = arg[2].replace("\n", "")

                table += "\n" + arg_template.format(
                    real_arg,
                    type_map.get(real_arg, "-"),
                    arg_desc,
                    rtype_map.get(real_arg, "-"),
                )
        fnsignature = self.docstring.split(":param")[0].strip()
        sig = self.signature.replace("\n", "").replace("\t", "").replace("    ", "")
        if arg_count:
            return f"### `{sig}`\n\n{fnsignature}\n{table}\n\n"
        return f"### `{sig}`\n\n{fnsignature}\n\n\n"


def extract_python_docs(file_text):
    for captured in c_joint_py.findall(file_text):
        print(captured)
        if captured:
            yield PythonDocstring(
                signature=captured[1].strip().replace("\n", ""),
                rtype=captured[2].strip().replace("\n", ""),
                docstring=captured[-1].strip(),
            )


def extract_cpp_docs(file_text):
    for captured in c_cdoc_rgx.findall(file_text):
        if captured:
            yield captured[1].strip()


def create_api_markdown_file(src_filename):
    _, file_extension = os.path.splitext(src_filename)
    target_filename = os.path.basename(os.path.dirname(src_filename)).replace(
        file_extension, ".md"
    )
    if not target_filename.endswith(".md"):
        target_filename += ".md"

    md_fn = ""
    with open(src_filename, "r") as f:
        ftext = f.read()

        class_docs = ftext.split("class")[1:]
        for i, doc_ in enumerate(class_docs):
            doc_ = (
                doc_.strip()
                .replace("@staticmethod", "")
                .replace("@classmethod", "")
                .replace("@overload", "")
            )
            class_name = doc_.partition("\n")[0].replace(":", "").strip()
            md_fn += f"## `{class_name}`"
            for g in extract_python_docs(doc_.replace("...", "...\n")):
                sig = g.py_signature_to_markdown()
                md_fn += f"\n{sig}\n"
            md_fn += "  \n"

    with open(
        os.path.join(os.path.dirname(__file__), GEN_FOLDER, target_filename), "w"
    ) as f:
        f.write(md_fn)


def on_startup(command, dirty, **kwargs):
    fn_lists = [
        *glob.glob(os.path.join(os.path.dirname(__file__), "..", "cmtj/*/*.pyi")),
        *glob.glob(os.path.join(os.path.dirname(__file__), "..", "cmtj/*.pyi")),
    ]
    for fn in fn_lists:
        create_api_markdown_file(fn)


if __name__ == "__main__":
    on_startup()
