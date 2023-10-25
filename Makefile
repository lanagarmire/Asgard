VERfile="VERSION.txt"
GBOXfile="GBOX_BASE_NAME.txt"
VER=`cat $(VERfile)`
GBOX=`cat $(GBOXfile)`:$(VER)
export

docker:
	docker build -t $(GBOX) .

docker-push:
	docker push $(GBOX)

server:
	docker run --rm -v `pwd`:/home/rstudio/Asgard -p 8787:8787 -it $(GBOX)

shell:
	docker run --rm -it $(GBOX) /bin/bash
