import streamlit as st
from PIL import Image

app_title = "CMTJ UI"
im = Image.open("assets/icon-larger.png")
st.logo(im)
domain_page = st.Page(
    "domain.py",
    title="Domain fitting",
    url_path="domains",
    icon=":material/bubble_chart:",
)
spectrum_page = st.Page(
    "streamlit_app.py",
    title="Spectrum fitting",
    url_path="spectrum",
    icon=":material/bar_chart:",
)

pg = st.navigation([spectrum_page, domain_page])
st.set_page_config(
    page_title=app_title,
    page_icon="🦈",
    layout="wide",
)
pg.run()
# > streamlit run ./view/streamlit_app.py --server.runOnSave true
