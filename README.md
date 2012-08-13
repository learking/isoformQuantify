# isoformQuantify: An interactive web service for isoform quantification using RNA-Seq data

isoformQuantify can be used to compute relative abundances of different transcripts 
belonging to the same gene and visualize the result.

## Installation of isoformQuantify

1. install Django (https://www.djangoproject.com/download/)
2. install scipy (http://www.scipy.org/SciPy)
3. install samtools (http://sourceforge.net/projects/samtools/files/samtools/0.1.18/)
4. clone the project repository to the location where you want to put
   
   git clone https://github.com/jaxcs/charlie.git

5. change line 20, 28, 29 in `/isoformQuantifyApp/getIsoformAbund.py` to appropriate directories in your server
6. copy data files in the `/isoformQuantifyApp/data/` directory to the same directory in the server (right now, data files have been ignored by git repo commit)
7. use command like `python manage.py runserver 0.0.0.0:8080` to kick off web service in the root folder that contains `manage.py`

## Directory structure

### Root directory
* `isoformQuantify`: configuration files for Django (Please refer to https://docs.djangoproject.com/en/1.4/)
* `isoformQuantifyApp`: actual project code locates

### `isoformQuantifyApp`
* `*.py`: Python scripts that implement algorithm and do scientific computation
* `/static/`: contains javascripts for visualization, css styling files
* `/templates/`: index.html (our homepage)

## Python scripts and their functions
* "views.py": respond to web request (which function responds to which request is determined by `/isoformQuantify/urls.py`)
* "models.py": stores Python objects used in algorithm
* "getIsoformAbund.py": calculate relative abundances of transcripts (main script that implements the algorithm)
* "test.py": auto generated file by Django (demonstrates how to write a test)

