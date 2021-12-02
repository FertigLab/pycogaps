FROM python:3.8

RUN useradd -ms /bin/bash user
USER user
WORKDIR /home/user

COPY requirements.txt .
COPY helper_functions.py .
COPY PyCoGAPS.py .
COPY setup.py .
COPY vignette.py .

RUN pip install -r requirements.txt

COPY src /src/ .

RUN python3 setup.py install
RUN python3 vignette.py

CMD ["python", "setup.py", "vignette.py"]


