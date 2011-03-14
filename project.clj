(defproject planar-quad "1.0.0-SNAPSHOT"
  :description "FIXME: write"
  :dependencies [[org.clojure/clojure "1.2.0"]
                 [org.clojure/clojure-contrib "1.2.0"]
		 [incanter/incanter-core "1.2.3-SNAPSHOT"]
                 [net.java.dev.jogl/jogl "1.1.1a"]]
  :native-dependencies [[net.java.dev.jogl/jogl-macosx-universal "1.1.1a"]
                        [org.clojars.angerman/jreality-native-deps "2011-01-29"]]
  :native-path "native/macosx/x86"
  :dev-dependencies [[native-deps "1.0.5"]
                     [swank-clojure "1.2.1"]])
