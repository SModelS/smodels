VERSION=`cat ../smodels/version`

all: tarballs

tarballs:
	./createTarballs.py

clean:
	rm -rf smodels-* 

version:
	@echo $(VERSION)

upload:
	scp ../ReleaseNotes smodels-v*.tgz smodels.hephy.at:/var/www/downloads/v1.1/
	scp -r smodels-v$(VERSION)/docs smodels.hephy.at:/var/www/docs/v$(VERSION)
