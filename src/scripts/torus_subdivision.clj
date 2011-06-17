(ns scripts.torus-subdivision
  (:require [jreality :as jr])
  (:require [geometry :as g])
  (:require [geometry.renderer :as r])
  (:require geometry.shapes.torus)
  (:use [geometry.subdivision.catmull-clark :only (catmull-clark-subd)])
  (:use [geometry.subdivision.doo-sabin :only (doo-sabin-subd)])
  (:use [geometry.subdivision.planar :only (planar-subd)])
  (:require [tikz :as t])
  (:require [clojure.contrib.duck-streams :as ds]))

(defn render-to-file [file geom]
  (let [f (t/resolve-file (t/temp-file file))
        eye  [5 6 4]
        up   [0 0 1]
        near 1.0
        far  13.0]
    (println (format "Writing file %s" f))
    (ds/with-out-writer f
      (println (t/header))
      (doseq [line (r/render-geometry geom eye up near far :magnify 15.0)]
        (println line))
      (println (t/footer)))))


(defn -main []
  (let [ifsf (jr/indexed-face-set-factory)
        geom (geometry.shapes.torus/create 10 5 1 3)]
    (render-to-file "torus.bare.tikz" geom)
    (g/with-geometry geom
      (planar-subd)
      ;;      (planar-subd)
      (render-to-file "torus.planar-subd.tikz" geom)
      (g/into-faceset-factory ifsf))
    (doall
     (let [geom (geometry.shapes.torus/create 10 5 1 3)]
       (g/with-geometry geom
         (doo-sabin-subd)
         (render-to-file "torus.doo-sabin-subd.tikz" geom))))
    (doall
     (let [geom (geometry.shapes.torus/create 10 5 1 3)]
       (g/with-geometry geom
         (catmull-clark-subd)
         (render-to-file "torus.catmull-clark-subd.tikz" geom))))
    (jr/show (jr/sgc ifsf))))
