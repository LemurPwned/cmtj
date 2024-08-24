import streamlit as st

app_title = "CMTJ UI"

domain_page = st.Page("domain.py", title="Domain fitting", url_path="domains", icon="📈")
spectrum_page = st.Page(
    "streamlit_app.py", title="Spectrum fitting", url_path="spectrum", icon="📈"
)

pg = st.navigation([domain_page, spectrum_page])

st.set_page_config(page_title=app_title, page_icon=":rocket:")
pg.run()
# > streamlit run ./view/streamlit_app.py --server.runOnSave true
