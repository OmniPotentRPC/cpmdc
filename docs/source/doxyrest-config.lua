INPUT_FILE = "xml/index.xml"
OUTPUT_FILE = "api/index.rst"
INDEX_TITLE = "API Reference"

local doxyrest_prefix = os.getenv("CONDA_PREFIX") or "/usr"

FRAME_DIR_LIST = {
    doxyrest_prefix .. "/share/doxyrest/frame/common",
    doxyrest_prefix .. "/share/doxyrest/frame/cfamily",
}

LANGUAGE = "cpp"
ESCAPE_ASTERISKS = true
ESCAPE_TRAILING_UNDERSCORES = true
