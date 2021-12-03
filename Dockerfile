FROM python:3.8

RUN useradd -ms /bin/bash user
USER user
WORKDIR /home/user

RUN git clone --recurse-submodules https://github.com/FertigLab/pycogaps.git

WORKDIR ./pycogaps/
RUN git pull

# COPY requirements.txt .
# COPY helper_functions.py .
# COPY subset_data.py .
# COPY PyCoGAPS.py .
# COPY setup.py .
# COPY vignette.py .

RUN pip install -r requirements.txt --user

# RUN aws s3 cp s3://pycogaps/GSE98638_HCC.TCell.S5063.count.txt ./data/
# RUN aws s3 cp s3://pycogaps/pheno.txt ./data/

# COPY src/ .
# COPY src/bindings.cpp ./src/
# COPY src/CoGAPS/ ./src/CoGAPS/

# COPY data/ .
# COPY data/GSE98638_HCC.TCell.S5063.count.txt ./data/
# COPY data/pheno.txt ./data/
# COPY data/GIST.csv ./data/

# WORKDIR /home/user

RUN python3 setup.py install --user
RUN python3 vignette.py

CMD ["python", "setup.py", "vignette.py"]


