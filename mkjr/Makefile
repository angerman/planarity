VERSION= $(shell date "+%Y-%m-%d")
DL=curl -s -O
EX=unzip -qq
build: clean pom jreality native cleanup
	@echo "Done..."

jreality:
	@echo "Downloading jReality libraries..."
	@rm -fR lib
	@mkdir lib
	@${DL} http://www3.math.tu-berlin.de/jreality/download/latest/jReality.jar
	@mv jReality.jar lib
	@${DL} http://www3.math.tu-berlin.de/jreality/download/latest/jReality-3rd-party.zip
	@${EX} jReality-3rd-party.zip

native:
	@echo "Downloading jReality Native libraries..."
	@rm -fR native
	@mkdir -p native/windows/x86
	@${DL} http://www3.math.tu-berlin.de/jreality/download/latest/jReality-jni-win32.zip
	@${EX} jReality-jni-win32.zip
	@mv jni/* native/windows/x86
	@mkdir -p native/windows/x86_64
	@${DL} http://www3.math.tu-berlin.de/jreality/download/latest/jReality-jni-win64.zip
	@${EX} jReality-jni-win64.zip
	@mv jni/* native/windows/x86_64
	@mkdir -p native/macosx/x86
	@${DL} http://www3.math.tu-berlin.de/jreality/download/latest/jReality-jni-macosx.zip
	@${EX} jReality-jni-macosx.zip
	@mv jni/* native/macosx/x86
	@mkdir -p native/linux/x86
	@${DL} http://www3.math.tu-berlin.de/jreality/download/latest/jReality-jni-linux32.zip
	@${EX} jReality-jni-linux32.zip
	@mv jni/* native/linux/x86
	@mkdir -p native/linux/x86_64
	@${DL} http://www3.math.tu-berlin.de/jreality/download/latest/jReality-jni-linux64.zip
	@${EX} jReality-jni-linux64.zip
	@mv jni/* native/linux/x86_64

cleanup:
	@echo "Cleanup..."
	@rm -fR jni
	@rm *.zip

clean:
	@rm -fR lib native pom.xml

pom:
	@rm -f pom.xml
	@echo "Creating POM file..."
	@echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" >> pom.xml
	@echo "<project>" >> pom.xml
	@echo "  <modelVersion>${VERSION}</modelVersion>" >> pom.xml
	@echo "  <groupId>org.clojars.${USER}</groupId>" >> pom.xml
	@echo "  <artifactId>jreality-native-deps</artifactId>" >> pom.xml
	@echo "  <version>${VERSION}</version>" >> pom.xml
	@echo "  <name>jReality</name>" >> pom.xml
	@echo " </project>" >> pom.xml

jar:
	@rm -f jreality-native-deps-${VERSION}.jar
	@echo "Creating jar..."
	@jar -cMf jreality-native-deps-${VERSION}.jar native lib

upload: jar
	@echo "Uploading to clojars..."
	@scp pom.xml jreality-native-deps-${VERSION}.jar clojars@clojars.org:
