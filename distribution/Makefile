VERSION=`cat ../smodels/version`

all: tarballs

tarballs:
	./createTarballs.py

clean:
	rm -rf smodels-*  create.log

version:
	@echo $(VERSION)

hidden:
	scp ../ReleaseNotes smodels*.tgz  smodels.hephy.at:/var/www/.hidden
	scp -r smodels-v$(VERSION)/docs smodels.hephy.at:/var/www/.hidden/docs/v$(VERSION)

upload:
	scp ../ReleaseNotes smodels*.tgz smodels.hephy.at:/var/www/downloads/v$(VERSION)/
	scp -r smodels-v$(VERSION)/docs/* smodels.hephy.at:/var/www/docs/v$(VERSION)
