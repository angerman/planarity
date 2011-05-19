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

(defn -main []
  (let [ifsf (jr/indexed-face-set-factory)
        geom (geometry.shapes.torus/create 10 5 1 3)]
    (g/with-geometry geom
      (planar-subd)
      (g/into-faceset-factory ifsf))
    (jr/show (jr/sgc ifsf))))
