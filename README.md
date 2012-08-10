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
6. use command like `python manage.py runserver 0.0.0.0:8080` to kick off web service in the root folder that contains `manage.py`

## 