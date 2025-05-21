"""Gradio interface for generating CMTJ simulation code with an LLM."""

from __future__ import annotations

import traceback

try:
    import gradio as gr
    from llama_index import VectorStoreIndex, SimpleDirectoryReader
    from llama_index.llms import OpenAI
    from memory_io import Memory
except Exception:  # pragma: no cover - used for offline placeholder
    # Provide minimal placeholders so the file can be imported without deps
    class Memory:
        def __init__(self):
            self.items = []

        def add(self, item: str) -> None:
            self.items.append(item)

        def history(self) -> str:
            return "\n".join(self.items)

    class Dummy:
        pass

    gr = Dummy()
    VectorStoreIndex = Dummy
    SimpleDirectoryReader = Dummy
    OpenAI = Dummy


ATTEMPTS = 3


def build_index() -> VectorStoreIndex:
    """Load repository docs and examples into LlamaIndex."""
    docs = []
    for path in ["docs", "examples", "python"]:
        try:
            reader = SimpleDirectoryReader(path)
            docs.extend(reader.load_data())
        except Exception:
            pass
    return VectorStoreIndex.from_documents(docs)


class CodeAgent:
    def __init__(self) -> None:
        self.memory = Memory()
        self.index = build_index()
        try:
            self.llm = OpenAI(model="gpt-3.5-turbo")
        except Exception:
            self.llm = None
        self.query_engine = self.index.as_query_engine()

    def __call__(self, prompt: str) -> str:
        errors: list[str] = []
        for _ in range(ATTEMPTS):
            context = ""
            try:
                context = self.query_engine.query(prompt).response
            except Exception:
                pass
            history = self.memory.history()
            sys_prompt = (
                f"Generate python code for the CMTJ library.\n"\
                f"Prompt: {prompt}\n"\
                f"Context: {context}\n"\
                f"Previous errors: {history}"
            )
            try:
                response = self.llm.complete(sys_prompt)
                code = response.text
            except Exception:
                # fallback if llm not available
                code = "# TODO: Generated code would appear here\n"
            try:
                exec(code, {})
                self.memory.add(code)
                return code
            except Exception as e:  # pragma: no cover - runtime execution check
                tb = traceback.format_exc()
                self.memory.add(tb)
                errors.append(tb)
        return f"Failed to generate working code after {ATTEMPTS} attempts. Last error: {errors[-1] if errors else 'n/a'}"


def main() -> None:
    agent = CodeAgent()
    iface = gr.Interface(fn=agent, inputs="text", outputs="code", title="CMTJ Simulation Code Generator")
    iface.launch()


if __name__ == "__main__":  # pragma: no cover - manual launch
    main()
