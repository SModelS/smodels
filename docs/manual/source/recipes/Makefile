all: convert clean # remove_ipynb

convert:
	./convert.py

clean:
	rm -f all_data_frame.txt # lheLLPExample.py

run:
	jupyter notebook
#ipython notebook

run_remotely:
	jupyter notebook --ip xx.xx.xx.xx --port 8888

remove_ipynb:
	rm -rf *.ipynb
