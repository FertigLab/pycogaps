FROM python:3.8

COPY . /pycogaps
WORKDIR pycogaps

RUN apt-get update && apt-get install -y libboost-all-dev
RUN pip install --upgrade pip 
RUN pip install -r requirements.txt

RUN python ./setup.py install

ENTRYPOINT ["python", "./vignette.py"]



