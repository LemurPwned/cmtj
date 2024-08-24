FROM python:3.11
EXPOSE 8501

WORKDIR /app
COPY view/ /app/
COPY assets/ /app/
RUN ls
RUN python3 -m pip install streamlit && \
    git clone https://github.com/LemurPwned/cmtj.git && \
    cd cmtj && \
    python3 -m pip install .[utils]

CMD streamlit run nav.py
