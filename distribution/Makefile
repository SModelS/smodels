VERSION=`cat ../smodels/version`

all: run

run:
	./createTarball.py

clean:
	rm -rf smodels-* 

version:
	@echo $(VERSION)

upload:
	scp ../ReleaseNotes smodels-v*.tgz smodels.hephy.at:/var/www/downloads/
	scp -r smodels-v$(VERSION)/docs smodels.hephy.at:/var/www/docs/v$(VERSION)
