all: run

run:
	./createTarball.py

clean:
	rm -rf smodels-v* 

upload:
	scp ../ReleaseNotes smodels-v*.tgz smodels.hephy.at:/var/www/downloads/
