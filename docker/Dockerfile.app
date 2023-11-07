FROM python:3.11
EXPOSE 8501

WORKDIR /app
COPY view/ /app/
RUN ls
RUN python3 -m pip install -r requirements.txt && \
    git clone https://github.com/LemurPwned/cmtj.git && \
    cd cmtj && \
    python3 -m pip install .[utils]

CMD streamlit run streamlit_app.py
