all: notebook

notebook:
	ipython notebook --matplotlib=inline

convert:
	ipython nbconvert --template full --to python example1.ipynb
	ipython nbconvert --template full --to html example1.ipynb
##ipython nbconvert --template basic --to html example1.ipynb
